#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#include "../RadialGrid/RadialGrid.hpp"

// a radial grid in lua is represented as a table that contains the radial grid pointer as light userdata
// and the appropriate metatable "LSMS.RadialGrid".
// radial grids that have been allocated by lua are considered temporary and contain the field "ownedByLua" set to true.
// only LSMS internal functions are allowed to assign pointer values and care should be taken that pointers owned by lua are never passed to such functions!

typedef struct
{
  RadialGrid *g;
  bool ownedByLua;
} RadialGridLua;

/*
static int newRadialGridLua(lua_State *L)
{
  RadialGridLua *t=(RadialGridLua *)lua_newuserdata(L,sizeof(RadialGridLua));
  luaL_getmetatable(L, "LSMS.RadialGrid");
  lua_setmetatable(L, -2);
  t->g=new RadialGrid;
  t->ownedByLua=true;
  return 1;
}
*/

static int newTemporaryRadialGridLua(lua_State *L)
{
  RadialGrid *t=new RadialGrid;
  lua_newtable(L);
  lua_pushlightuserdata(L,t);
  lua_setfield(L,-2,"RadialGrid");
  lua_pushboolean(L,1);
  lua_setfield(L,-2,"ownedByLua");

  luaL_getmetatable(L, "LSMS.RadialGrid");
  lua_setmetatable(L, -2);

  return 1;
}

/*
static int copyRadialGridLua(lua_State *L)
{
  RadialGridLua *g=(RadialGridLua *)lua_touserdata(L,1);
  RadialGridLua *t=(RadialGridLua *)lua_newuserdata(L,sizeof(RadialGridLua));
  luaL_getmetatable(L, "LSMS.RadialGrid");
  lua_setmetatable(L, -2);
  t->g=new RadialGrid;
  *(t->g) = *(g->g);
  t->ownedByLua=true;
  return 1;
}
*/

static int radialGridLua__gc(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  lua_getfield(L,1,"ownedByLua");
  if(lua_toboolean(L,-1)) delete g;
  return 0;
}

static int generateRadialGridLua(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  Real x0=luaL_checknumber(L,2);
  Real h=luaL_checknumber(L,3);
  int N=luaL_checkint(L,4);
  int jmt=luaL_checkint(L,5);
  int jws=luaL_checkint(L,6);

  generateRadialGrid(g,x0,h,N,jmt,jws);

  return 0;
}

static int sizeRadialGridLua(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  lua_pushnumber(L, g->N);

  return 1;
}

static int getRadialGridLua(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  int idx = luaL_checkint(L,2);

  luaL_argcheck(L, 0<=idx && idx < g->N,2,"index out of range");
  lua_pushnumber(L,g->r_mesh[idx]);
  return 1;
}

static int getRadialGridXLua(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  int idx = luaL_checkint(L,2);

  luaL_argcheck(L, 0<=idx && idx < g->N,2,"index out of range");
  lua_pushnumber(L,g->x_mesh[idx]);
  return 1;
}

static int luaRadialGridToString(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  if(g->N > 0)
    lua_pushfstring(L, "RadialGrid(%f,%f,%d,%d,%d)", g->x_mesh[0],
                    g->h, g->N, g->jmt, g->jws);
  else
    lua_pushfstring(L, "RadialGrid(0)");
  return 1;
}

static int radialGridLua__index(lua_State *L)
{
  lua_getfield(L,1,"RadialGrid");
  RadialGrid *g=(RadialGrid *)lua_touserdata(L,-1);
  if(lua_type(L,2)==LUA_TNUMBER)
    return getRadialGridLua(L);
  lua_getmetatable(L,1);
  lua_pushvalue(L,2);
  lua_gettable(L,-2);
  return 1;
}

static const struct luaL_Reg RadialGrid_lib_f [] = {
  {"new", newTemporaryRadialGridLua},
  {NULL, NULL}
};

static const struct luaL_Reg RadialGrid_lib_m [] = {
  {"__index", radialGridLua__index},
  {"generate", generateRadialGridLua},
  {"__tostring", luaRadialGridToString},
  {"__gc", radialGridLua__gc},
  {"size", sizeRadialGridLua},
  {"get", getRadialGridLua},
  {"x", getRadialGridXLua},
//  {"copy", copyRadialGridLua},
  {NULL, NULL}
};

int luaopen_RadialGrid (lua_State *L)
{
  luaL_newmetatable(L, "LSMS.RadialGrid");
  luaL_register(L, NULL, RadialGrid_lib_m);

  // lua_pushstring(L, "__index");
  // lua_pushstring(L,"get");
  // lua_gettable(L,2);
  // lua_settable(L,1);

  luaL_register(L, "RadialGrid", RadialGrid_lib_f);
  return 1;
}
