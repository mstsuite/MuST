#ifndef LSMS_LUASUPPORT_H
#define LSMS_LUASUPPORT_H

#include "lua.hpp"
//#include "lua.h"
#include "Real.hpp"

void luaStackDump(lua_State *L);

bool luaGetStrN(lua_State *L, const char *name, char *s, int n);
bool luaGetReal(lua_State *L, const char *name, Real *val);
bool luaGetInteger(lua_State *L, const char *name, int *val);
bool luaGetBoolean(lua_State *L, const char *name, bool *val);

bool luaGetFieldInTable(lua_State *L, const char *name, const char *field);
bool luaGetFieldFromStack(lua_State *L, const char *field);
bool luaGetPositionInTable(lua_State *L, const char *name, int idx);
bool luaGetRealFieldInTable(lua_State *L, const char *name, const char *field, Real *val);
bool luaGetIntegerFieldInTable(lua_State *L, const char *name, const char *field, int *val);
bool luaGetRealPositionInTable(lua_State *L, const char *name, int idx, Real *val);
bool luaGetIntegerPositionInTable(lua_State *L, const char *name, int idx, int *val);
bool luaGetStrNFromStack(lua_State *L, const char *name, char *s, int n);
bool luaGetRealPositionFromStack(lua_State *L, int idx, Real *val);
bool luaGetRealFieldFromStack(lua_State *L, const char *field, Real *val);
bool luaGetIntegerFieldFromStack(lua_State *L, const char *field, int *val);
#endif
