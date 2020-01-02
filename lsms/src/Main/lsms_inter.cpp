#include <stdio.h>
#include <string.h>
#include "lua.hpp"
//#include "lua.h"
//#include "lauxlib.h"
//#include "lualib.h"

void initLSMSLuaInterface(lua_State *L);

int main(int argc, char *argv[])
{
  char buff[256];
  int error;

  lua_State *L=lua_open();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  while(fgets(buff,sizeof(buff), stdin) != NULL)
  {
    error = luaL_loadbuffer(L, buff, strlen(buff), "line") ||
            lua_pcall(L, 0, 0, 0);
    if(error)
    {
      fprintf(stderr, "%s", lua_tostring(L, -1));
      lua_pop(L,1);
    }
  }

  lua_close(L);
  return 0;
}
