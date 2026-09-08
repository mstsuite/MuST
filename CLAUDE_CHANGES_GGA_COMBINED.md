# `MuST_new_GGA` — combined change set relative to `MuST`

Companion to `CLAUDE_CHANGES_GGA_COMBINED.patch` in this directory. That patch is
the *only* content difference between this tree and `~/Desktop/MuST`, and this
file is the *only* other file added.

---

## 1. Provenance and how the two inputs were combined

| Input | Commits | Base |
|---|---|---|
| `CLAUDE_CHANGES.patch` | 1 | `2617706e` |
| `CLAUDE_CHANGES_GGA.patch` | 3 (`1/3`, `2/3`, `3/3`) | `2617706e` |

**The two inputs are not independent.** Commit `1/3` of the GGA series is
byte-identical to the single commit of `CLAUDE_CHANGES.patch`; the only textual
differences are the `format-patch` subject counter (`[PATCH]` vs `[PATCH 1/3]`)
and one trailing newline:

```
$ diff <(sed -n '1,1442p' CLAUDE_CHANGES_GGA.patch) CLAUDE_CHANGES.patch
4c4
< Subject: [PATCH 1/3] GPU-accelerate the L-expanded density interpolation in
---
> Subject: [PATCH] GPU-accelerate the L-expanded density interpolation in
1441d1440
<
```

So `CLAUDE_CHANGES_GGA` **supersedes** `CLAUDE_CHANGES` rather than
complementing it. Combining is therefore a *squash of the 3-commit GGA series*,
not a merge of two independent change sets. Applying `CLAUDE_CHANGES.patch` in
addition would be a no-op at best and a conflict at worst.

The three squashed commits:

1. `GPU-accelerate the L-expanded density interpolation in calExchangeJl`
   — value-only batching (the `B1` hotspot).
2. `Link DensityOnGridModule into the KUBO executable`
   — one line in `KUBO/src/Makefile`; without it `KUBO` fails to link once
   `PotentialGenerationModule` references the new module.
3. `GGA: batch the density gradient on the same GPU/DGEMM path`
   — extends the batching to $\nabla\rho$, which set (1) explicitly left on the
   per-point path.

### Construction and equivalence check

The series was split with `csplit` on the `format-patch` commit boundary and
applied sequentially with `patch -p1` to a pristine copy of the eight touched
files; a single squashed diff was then generated and re-applied to a second
pristine copy. The two results are identical:

```
$ diff -r sequential/ squashed/     # no output
```

Each of the eight files in this tree was verified by `md5sum` against the
sequentially-patched reference. `patch` reported no fuzz, no offsets and no
`.rej`/`.orig` files.

### Base-drift check

The patches were authored against `2617706e`; `MuST` is now at `69b1342b`
(three commits later: `824a4778`, `e7ede034`, `69b1342b`). All six *modified*
files are unchanged across that interval —

```
$ git diff --stat 2617706e..HEAD -- <each of the six>   # empty for all six
```

— so the series applies to the current `MuST` without adaptation, and the two
*new* files do not collide with anything added since the base.

---

## 2. Diffstat

```
 KUBO/src/Makefile                      |   1 +
 MST/Accelerator/CMakeLists.txt         |   5 +-
 MST/Accelerator/DensityInterp_Accel.cu | 883 +++++++++++++++++++++++++++++++++
 MST/Accelerator/Makefile               |   1 +
 MST/src/AngularIntegrationModule.F90   |  53 +-
 MST/src/DensityOnGridModule.F90        | 785 +++++++++++++++++++++++++++++
 MST/src/Makefile                       |   1 +
 MST/src/PotentialGenerationModule.F90  | 155 +++++-
 8 files changed, 1877 insertions(+), 7 deletions(-)
```

Two new files, six modified. No file deleted, no file renamed.

---

## 3. What the change set does

### 3.1 The cost being removed

`PotentialGenerationModule::calExchangeJl` evaluated the density pointwise via
`ChargeDensityModule::getChargeDensityAtPoint`. Each call performed a `hunt()`
bisection over the radial mesh plus an `n_inter`-point Neville interpolation for
every $jl$ component. The targets, however, are

$$\mathbf{r}_{i,g} = r_i\,\hat{u}_g ,\qquad r_i \in \{r\text{-mesh nodes}\},$$

i.e. exactly on mesh nodes, and the directions $\hat u_g$ are fixed for the whole
run (`AngularIntegrationModule::setAngularData` builds them once). With
$n_r \sim 1500$ and $n_g = 50\times 80 = 4000$ that is $\sim 6\times10^6$ point
evaluations per (atom, species) per SCF iteration.

### 3.2 Factorization — value

Radial interpolation is **linear** in the stored coefficients, so it separates
from the angular sum:

$$A(i,jl) \;=\; \sum_{k=1}^{n_\mathrm{inter}} w_k(r_i)\,\rho_L\!\big(\mathrm{irp}(i)+k-1,\,jl\big)$$

$$\rho(r_i,\hat u_g) \;=\; \sum_{jl} f_{a2}(jl)\,\mathrm{Re}\!\left[A(i,jl)\,Y_{jl}(\hat u_g)\right] \;=\; \big[A_r W_r - A_i W_i\big](i,g)$$

with

$$W_r(jl,g) = f_{a2}(jl)\,\mathrm{Re}\,Y_{jl}(\hat u_g),\qquad W_i(jl,g) = f_{a2}(jl)\,\mathrm{Im}\,Y_{jl}(\hat u_g),$$

built **once**. The angular sum is then two real GEMMs.

The interpolation weights $w_k$ are **Lagrange** weights, not the Neville
recursion — same interpolating polynomial, but linear in the data, which is what
permits the factorization above.

### 3.3 Factorization — gradient (the GGA increment)

In `SphericalHarmonicsModule::SphericalHarmonics4` the gradient scaling is
`rfac = clm(jl)/r`, and the $\hat\theta,\hat\varphi$ unit vectors depend on
direction only. The radius therefore enters `grad_ylm` **only** as an overall
$1/r$, so the angular part can be tabulated once on the unit sphere, exactly like
$Y_{jl}$. With $\dot A(i,jl)$ the interpolated radial derivatives:

$$D(i,g) \;=\; \sum_{jl} f_{a2}(jl)\,\mathrm{Re}\!\left[\dot A(i,jl)\,Y_{jl}(\hat u_g)\right]$$

$$T_c(i,g) \;=\; \sum_{jl} f_{a2}(jl)\,\mathrm{Re}\!\left[A(i,jl)\,G_{jl,c}(\hat u_g)\right],\qquad G_{jl,c} = f_{a2}(jl)\left[r\,\partial_{x_c} Y_{jl}\right]_{r=1}$$

$$\boxed{\;\partial_c\rho(i,g) \;=\; D(i,g)\,\hat u_c(g) \;+\; \frac{T_c(i,g)}{r_i}\;}$$

This reproduces term by term the CPU expression in `getChargeDensityAtPoint`,

```
grad(i) += fa2(jl)*Re[ der_rho_in*ylm(kl)*er(i) + rho_in(jl)*grad_ylm(kl,i) ]
```

with $\hat e_r = \hat u_g$.

$D$ reuses the existing $W_r,W_i$ tables (one extra GEMM pair); $T_c$ needs three
new tables `Gr(:,:,c)`, `Gi(:,:,c)` and three more GEMM pairs. Total **eight**
GEMM pairs per (atom, species) instead of two, against
$\mathcal{O}(n_r\,n_g\,j_{\max})$ Neville interpolations before.

---

## 4. Symbol → code map

| Symbol | Code | Location |
|---|---|---|
| $\rho_L(ir,jl)$ | `p_den_l`, from `getChargeDensity("TotalNew",id,ia)` | `PotentialGenerationModule` |
| $\partial_r\rho_L(ir,jl)$ | `p_der_den_l`, 4th arg of `getChargeDensity` | `PotentialGenerationModule` |
| $A(i,jl)$ | `Ar`, `Ai` / `d_Ar`, `d_Ai` | `DensityOnGridModule`, `DensityInterp_Accel.cu` |
| $\dot A(i,jl)$ | same buffers, derivative pass | `DensityInterp_Accel.cu` |
| $W_r,W_i$ | `Wr`, `Wi` / `d_Wr`, `d_Wi` | built in `initDensityOnGrid` |
| $G_{jl,c}$ | `Gr(:,:,c)`, `Gi(:,:,c)` / `d_Gr[c]`, `d_Gi[c]` | built when `needGrad=.true.` |
| $D(i,g)$ | `d_D` | `DensityInterp_Accel.cu` |
| $T_c(i,g)$ | `d_T` (one buffer, reused for `c=0,1,2`) | `DensityInterp_Accel.cu` |
| $\rho(r_i,\hat u_g)$ | `den_grid(ir,ing)` | `calExchangeJl` |
| $m(r_i,\hat u_g)$ | `mom_grid(ir,ing)` | `calExchangeJl` |
| $\partial_c\rho$ | `dgrad_grid(ir,(c-1)*ngl+ing)` | `calExchangeJl` |
| $\partial_c m$ | `mgrad_grid(ir,(c-1)*ngl+ing)` | `calExchangeJl` |
| $w_k(r_i)$ | `wgt_tab`, Lagrange weights | `DensityOnGridModule` |
| $\mathrm{irp}(i)$ | `irp_tab(i)` | `DensityOnGridModule` |
| $\hat u_g$ | `upos(1:3,ing)` / `d_upos` | `AngularIntegrationModule` |

---

## 5. File-by-file delta

### New: `MST/Accelerator/DensityInterp_Accel.cu` (883 lines)

Eight `extern "C"` entry points, trailing-underscore / by-pointer convention of
`LSMS_Accel.cu`, own stream (`di_stream`) and cuBLAS handle, rank-modulo-device
assignment as in `init_lsms_gpu`:

| Entry point | Args | Role |
|---|---|---|
| `init_density_interp_gpu_` | 7 | allocate device buffers, create stream/handle |
| `push_angular_ylm_gpu_` | 6 | pack and upload $W_r,W_i$ |
| `push_angular_gradylm_gpu_` | 7 | pack and upload $G_r,G_i$ and $\hat u$ |
| `push_radial_mesh_gpu_` | 2 | upload radial mesh (cached by `mesh_id`) |
| `push_density_l_gpu_` | 5 | upload $\rho_L$ / $\partial_r\rho_L$ into a field slot |
| `eval_density_sphere_gpu_` | 5 | value only |
| `eval_density_grad_sphere_gpu_` | 7 | value and gradient in one call |
| `finalize_density_interp_gpu_` | 0 | free everything |

Kernels: device bracket search + Lagrange weights, strided $L$-component
interpolation, `di_packYlmWeightsKernel`, `di_packGradYlmWeightsKernel`, and
`di_assembleGradKernel` ($T_c \leftarrow D\hat u_c + T_c/r$, in place).

Ordering inside `eval_density_grad_sphere_gpu_` is load-bearing: the
**derivative** pass runs first and parks its contraction in `d_D`, freeing
`d_Ar`/`d_Ai` for the value pass. The reverse order would require a
device→host→device round trip of $n_r n_g$ doubles per component. `d_T` is a
single buffer reused for $c=0,1,2$; correctness rests on every kernel, cuBLAS
call and 2-D copy being issued on the single stream `di_stream`, which
serializes them.

### New: `MST/src/DensityOnGridModule.F90` (785 lines)

Public interface:

```fortran
public :: initDensityOnGrid, endDensityOnGrid,                  &
          calDensityOnAngularGrid, calDensityGradOnAngularGrid,  &
          isDensityOnGridGPU, isDensityGradAvailable,            &
          getDensityOnGridTime, printDensityOnGridInfo,          &
          FIELD_CHARGE, FIELD_MOMENT,                            &
          FIELD_DER_CHARGE, FIELD_DER_MOMENT
```

`FIELD_CHARGE = 1`, `FIELD_MOMENT = 2`, `FIELD_DER_CHARGE = 3`,
`FIELD_DER_MOMENT = 4` — the last two are the GGA increment (device buffer slots
3 and 4). `initDensityOnGrid(nrmax, jmaxmax, lmax, mype, iprint, needGrad)`
builds the gradient tables only when `needGrad = .true.`.

The CPU path uses BLAS `DGEMM` and is **algebraically identical** to the GPU
path, so the restructuring pays off with `ACCEL` off and the two paths can be
cross-checked against each other.

The unit-sphere tables come from `calYlm(upos(:,g), lmax, ylm_u, grady_u)` where
`upos` is `AngularIntegrationModule`'s own unit-vector table, so the directions
are bit-identical to those the point loop obtains from `getUnitVec`.

### Modified: `MST/src/PotentialGenerationModule.F90` (+152 / −3)

Four hunks:

* `initPotentialGeneration` — `initDensityOnGrid(jend_max, jmax_max, lmax_max, MyPEinGroup, node_print_level, needGrad=gga_functional)`, placed **after** `initAngularIntegration` because it caches that module's Ylm table.
* `endPotentialGeneration` — `endDensityOnGrid()` before `endAngularIntegration()`.
* `calExchangeJl` — batch evaluation hoisted above the `ing`/`ir` loops, gated by

  ```fortran
  use_batched = (.not.gga_functional) .or. isDensityGradAvailable()
  ```

  For GGA the call is `calDensityGradOnAngularGrid`, filling
  `dgrad_grid(jend, ngl*3)` (and `mgrad_grid` when `n_spin_pola == 2`) with
  component $c$ in columns `(c-1)*ngl+1 : c*ngl`. The point loop then reads
  `den_grid`/`dgrad_grid` instead of calling `getChargeDensityAtPoint`, and
  dispatches to `calExchangeCorrelation` exactly as before.
* deallocation / `nullify` after the loops.

The per-point fallback is retained verbatim on the `else` branches, so no
configuration silently loses its gradient.

### Modified: `MST/src/AngularIntegrationModule.F90` (+51 / −2)

Adds `getYlmTable()` and `getUnitVecTable()` to the `public` list and gives
`ylm` and `upos` the `target` attribute so pointers can be returned. Both accessors
guard on `Initialized` via `ErrorHandler`. No change to existing behaviour.

### Modified: build files (+1 line each, +3/−2 for CMake)

| File | Change |
|---|---|
| `MST/Accelerator/Makefile` | `DensityInterp_Accel.o` added to `OBJ0` |
| `MST/Accelerator/CMakeLists.txt` | `DensityInterp_Accel.cu` added to the source list |
| `MST/src/Makefile` | `DensityOnGridModule.o` added |
| `KUBO/src/Makefile` | `$(MST_ODIR)/DensityOnGridModule.o` added |

---

## 6. Where the result is exactly the old result

For on-node targets one Lagrange numerator factor is exactly `0.0` and the
surviving ratio is a quotient of two identical products, hence exactly `1.0`. The
interpolation is therefore **bitwise identical** to the Neville path at the mesh
nodes — which is where every target lies. The remaining differences are the
summation order of the angular contraction (GEMM vs scalar loop), so agreement
with the previous code should be at round-off, not merely at tolerance.

---

## 7. Cost

| | Before | After |
|---|---|---|
| Radial work per (atom, species) | $\mathcal{O}(n_r n_g j_{\max})$ Neville interps | $\mathcal{O}(n_r n_\mathrm{inter} j_{\max})$ |
| Angular work | fused into the above | 2 GEMM pairs (LDA), 8 GEMM pairs (GGA) |
| Extra host memory | — | $W_r,W_i$: $j_{\max}\times n_g$; $G_r,G_i$: $3\,j_{\max}\times n_g$ |
| Extra device memory (GGA) | — | `d_D`, `d_T`: $2\times n_r n_g$ doubles ≈ 48 MB each at $n_r=1500$, $n_g=4000$ |

Gradient tables and buffers are allocated only when `needGrad = .true.`.

---

## 8. Not verified — read before running

1. **Nothing here is compile-verified.** Neither the original authoring
   environment nor this one has `gfortran` or `nvcc`. Verification performed was
   structural only, and is reported here as exactly that:
   * brace / paren / bracket balance in `DensityInterp_Accel.cu` — 78/78, 459/459, 64/64;
   * `do`/`enddo`, `if…then`/`endif`, `subroutine`/`end subroutine` pairing in the three
     `.F90` files, each compared against the *unpatched* file so that
     regex artefacts cancel — no imbalance introduced;
   * every Fortran call site cross-checked against the CUDA entry-point arity —
     13 call sites, all matching (`init` 7/7, `push_angular_ylm` 6/6,
     `push_angular_gradylm` 7/7, `push_radial_mesh` 2/2 ×4, `push_density_l` 5/5 ×3,
     `eval_density_sphere` 5/5, `eval_density_grad_sphere` 7/7, `finalize` 0/0);
   * `use … only :` imports in `PotentialGenerationModule` all present in
     `DensityOnGridModule`'s `public` list.

   None of this substitutes for compiling. **Compile before running.**

2. **Open question — radial range when `size(p_den_l,1) < jend`.**
   `calExchangeJl` sets

   ```fortran
   nr_den = min(jend, size(p_den_l,1))
   ```

   and `calDensityOnAngularGrid` writes `den_grid(i,g)` only for `i = 1..nr`.
   Rows `nr_den+1 … jend` therefore stay at their initialized `ZERO`, and the
   point loop's `if (rho <= ZERO) cycle Loop_ir` skips them. On the old path
   those radial points went through `getChargeDensityAtPoint`, which clamps its
   stencil and returns a value. **If `size(p_den_l,1)` can be smaller than
   `Grid%jend` for any atom, this silently drops radial points.** Unverified in
   either direction — check whether `getChargeDensity("TotalNew",…)` is always
   dimensioned to at least `Grid%jend`.

3. **Numerical cross-check to run first.** A non-GGA full-potential run should
   reproduce the previous `PotNL_*` output and total energy to round-off. For
   GGA, the reference is the *value-only* tree (`MuST` + commit `1/3` only),
   which computes $\nabla\rho$ by the original per-point route; the two should
   agree very tightly since interpolation on mesh nodes is exact. Also worth
   running: `ACCEL` on vs off in this same tree, which exercises the GPU and CPU
   paths against each other.

4. **`nr < n_inter` and `nr > nr_max` are hard `ErrorHandler` stops** inside
   `calDensityOnAngularGrid`, as are `ld_grid < nr` and `ld_den < nr`. An atom
   with fewer than `n_inter` radial points would now abort where it previously
   ran.

---

## 9. Build notes

* The GPU path compiles in only when the accelerator is enabled: the top-level
  `Makefile` adds `-DACCEL -DCUDA` to `FPPFLAGS` when `Acceleration = 1`. See
  `arch/bolt_intel_accel` or `arch/summit_pgi_accel` for a working
  `ACCEL_CXX` / `ACCEL_OPT` / `ADD_LIBS` set — `-lcudart -lcuda -lcublas` are
  already linked there.
* With `Acceleration` unset, `DensityInterp_Accel.o` is `touch`-ed empty by the
  existing `Accelerator/Makefile` rule and `DensityOnGridModule` takes its CPU
  DGEMM path. No CUDA toolchain required.

---

## 10. Scope — what is *not* here

* Exchange-correlation evaluation itself remains on the CPU, single threaded.
  OpenMP there requires a thread-safety audit of `ExchCorrFunctionalModule`
  first: `calExchangeCorrelation` and `getExchCorrPot` / `getExchCorrEnDen`
  communicate through module state, so the loop cannot be wrapped in
  `!$omp parallel do` as written.
* `MST_Improvement_Notes` items **A1** (`calIntraPot` `nRpts+1` / `nRpts_ps+1`
  size mismatch), **A2** (last-channel-wins component flag in
  `computeNewPotential`), **A3** (`cycle` vs `exit` in `setPotComponentFlag`),
  **B2** (`calInterPlusMadPot`), **B3** (`calPseudoDipoleField` loop order),
  **B4** (OpenMP) and **B5** (dead `XchgCorrHat*` storage) are untouched.

---

## 11. Tree provenance

`MuST_new_GGA` is a copy of `~/Desktop/MuST` at `HEAD = 69b1342b`, working tree
clean, with `.git/` excluded (326 MB of history not duplicated; no commits,
branches or pushes were made in either tree). Verified with

```
$ rsync -ain --delete --exclude='.git/' MuST/ MuST_new_GGA/
```

which, before the patch was applied, reported no differences over all 15 322
files. After the patch the tree contains 15 325 files: the 15 322 originals with
six modified, plus `DensityInterp_Accel.cu`, `DensityOnGridModule.F90`, and the
two files described here.
