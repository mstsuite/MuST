#include <string.h>
#include "Real.hpp"
#include "LuaSupport.hpp"

// useful routine borrowed from internet
void luaStackDump(lua_State *L) {

  printf("Lua Stack: ");
  int i;
  int top = lua_gettop(L);
  for (i = 1; i <= top; i++) {  /* repeat for each level */
    int t = lua_type(L, i);
    switch (t) {

      case LUA_TSTRING:  /* strings */
        printf("`%s'", lua_tostring(L, i));
        break;

      case LUA_TBOOLEAN:  /* booleans */
        printf(lua_toboolean(L, i) ? "true" : "false");
        break;

      case LUA_TNUMBER:  /* numbers */
        printf("%g", lua_tonumber(L, i));
        break;

      default:  /* other values */
        printf("%s", lua_typename(L, t));
        break;

    }
    printf("  ");  /* put a separator */
  }
  printf("\n");  /* end the listing */
  fflush(stdout);
}

bool luaGetStrN(lua_State *L, const char *name, char *s, int n)
{
  lua_getglobal(L,name);
  if(!lua_isstring(L,-1)) {lua_pop(L,1); return false;}
  strncpy(s, lua_tostring(L,-1), n);
  lua_pop(L,1);
  return true;
}

bool luaGetReal(lua_State *L, const char *name, Real *val)
{
  lua_getglobal(L,name);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetInteger(lua_State *L, const char *name, int *val)
{
  lua_getglobal(L,name);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetBoolean(lua_State *L, const char *name, bool *val)
{
  lua_getglobal(L,name);
  if(!lua_isboolean(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_toboolean(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetFieldInTable(lua_State *L, const char *name, const char *field)
{
  lua_getglobal(L,name);
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  return true;
}

bool luaGetFieldFromStack(lua_State *L, const char *field)
{
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  return true;
}

bool luaGetPositionInTable(lua_State *L, const char *name, int idx)
{
  lua_getglobal(L,name);
  if(lua_isnil(L,-1)) {lua_pop(L,1); return false;}
  lua_pushinteger(L,idx);
  lua_gettable(L,-2);
  return true;
}

bool luaGetRealFieldInTable(lua_State *L, const char *name, const char *field, Real *val)
{
  if(!luaGetFieldInTable(L,name,field)) return false;
  if(!lua_isnumber(L,-1)) {lua_pop(L,2); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,2);
  return true;
}

bool luaGetIntegerFieldInTable(lua_State *L, const char *name, const char *field, int *val)
{
  if(!luaGetFieldInTable(L,name,field)) return false;
  if(!lua_isnumber(L,-1)) {lua_pop(L,2); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,2);
  return true;
}

bool luaGetRealPositionInTable(lua_State *L, const char *name, int idx, Real *val)
{
  if(!luaGetPositionInTable(L,name,idx)) return false;
  if(!lua_isnumber(L,-1)) {lua_pop(L,2); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,2);
  return true;
}

bool luaGetIntegerPositionInTable(lua_State *L, const char *name, int idx, int *val)
{
  if(!luaGetPositionInTable(L,name,idx)) return false;
  if(!lua_isnumber(L,-1)) {lua_pop(L,2); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,2);
  return true;
}

bool luaGetRealPositionFromStack(lua_State *L, int idx, Real *val)
{
  lua_pushinteger(L,idx);
  lua_gettable(L,-2);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetStrNFromStack(lua_State *L, const char *name, char *s, int n)
{
  lua_pushstring(L,name);
  lua_gettable(L,-2);
  if(!lua_isstring(L,-1)) {lua_pop(L,1); return false;}
  strncpy(s, lua_tostring(L,-1), n);
  lua_pop(L,1);
  return true;
}

bool luaGetRealFieldFromStack(lua_State *L, const char *field, Real *val)
{
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetIntegerFieldFromStack(lua_State *L, const char *field, int *val)
{
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,1);
  return true;
}
