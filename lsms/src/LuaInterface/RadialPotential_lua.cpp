#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#include "../RadialGrid/RadialGrid.hpp"
#include "../RadialGrid/RadialPotential.hpp"

// a radial grid in lua is represented as a table that contains the radial grid pointer as light userdata
// and the appropriate metatable "LSMS.RadialGrid".
// radial grids that have been allocated by lua are considered temporary and contain the field "ownedByLua" set to true.
// only LSMS internal functions are allowed to assign pointer values and care should be taken that pointers owned by lua are never passed to such functions!


static int newTemporaryRadialPotentialLua(lua_State *L)
{
  RadialPotential *t=new RadialPotential;
  lua_newtable(L);
  lua_pushlightuserdata(L,t);
  lua_setfield(L,-2,"RadialPotential");
  lua_pushboolean(L,1);
  lua_setfield(L,-2,"ownedByLua");

  luaL_getmetatable(L, "LSMS.RadialPotential");
  lua_setmetatable(L, -2);

  return 1;
}

static int radialPotentialLua__gc(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *g=(RadialPotential *)lua_touserdata(L,-1);
  lua_getfield(L,1,"ownedByLua");
  if(lua_toboolean(L,-1)) delete g;
  return 0;
}

static int getRadialGridFromPotentialLua(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *p=(RadialPotential *)lua_touserdata(L,-1);
  RadialGrid *g = p->g;

  lua_newtable(L);
  lua_pushlightuserdata(L,g);
  lua_setfield(L,-2,"RadialGrid");
//  lua_pushboolean(L,1);
//  lua_setfield(L,-2,"ownedByLua");

  luaL_getmetatable(L, "LSMS.RadialGrid");
  lua_setmetatable(L, -2);

  return 1;
}


static int getRadialPotentialLua(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *p=(RadialPotential *)lua_touserdata(L,-1);
  int is = luaL_checkint(L,2);
  int ir = luaL_checkint(L,3);

  luaL_argcheck(L, 0<=ir && ir < p->g->N,2,"index out of range");
  lua_pushnumber(L,p->vr(ir,is));
  return 1;
}

static int syncRadialPotentialLua(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *p=(RadialPotential *)lua_touserdata(L,-1);
  p->sync();
  return 0;
}

static int sizeRadialPotentialLua(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *g=(RadialPotential *)lua_touserdata(L,-1);
  lua_pushnumber(L, g->vr.n_row());

  return 1;
}

static int luaRadialPotentialToString(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *g=(RadialPotential *)lua_touserdata(L,-1);
//  if(g->N > 0)
//    lua_pushfstring(L, "RadialPotential(%f,%f,%d,%d,%d)", g->x_mesh[0],
//                    g->h, g->N, g->jmt, g->jws);
//  else
  lua_pushfstring(L, "RadialPotential(%p)",g);
  return 1;
}

static int radialPotentialLua__index(lua_State *L)
{
  lua_getfield(L,1,"RadialPotential");
  RadialPotential *g=(RadialPotential *)lua_touserdata(L,-1);
  //if(lua_type(L,2)==LUA_TNUMBER)
  //  return getRadialGridLua(L);
  lua_getmetatable(L,1);
  lua_pushvalue(L,2);
  lua_gettable(L,-2);
  return 1;
}

static const struct luaL_Reg RadialPotential_lib_f [] = {
  {"new", newTemporaryRadialPotentialLua},
  {NULL, NULL}
};

static const struct luaL_Reg RadialPotential_lib_m [] = {
  {"__index", radialPotentialLua__index},
//  {"generate", generateRadialGridLua},
  {"__tostring", luaRadialPotentialToString},
  {"__gc", radialPotentialLua__gc},
  {"grid", getRadialGridFromPotentialLua},
  {"size", sizeRadialPotentialLua},
  {"get", getRadialPotentialLua},
  {"sync", syncRadialPotentialLua},
//  {"copy", copyRadialGridLua},
  {NULL, NULL}
};

int luaopen_RadialPotential (lua_State *L)
{
  luaL_newmetatable(L, "LSMS.RadialPotential");
  luaL_register(L, NULL, RadialPotential_lib_m);

  // lua_pushstring(L, "__index");
  // lua_pushstring(L,"get");
  // lua_gettable(L,2);
  // lua_settable(L,1);

  luaL_register(L, "RadialPotential", RadialPotential_lib_f);
  return 1;
}
