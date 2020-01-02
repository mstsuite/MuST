/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_LSMSMODE_H
#define LSMS_LSMSMODE_H

#define MY_LSMSMODE_ENUM \
X(main) \
X(liz0)

#define X(name) name,
 
#define MY_LSMSMODE_ENUM_NAME LSMSMode

enum class MY_LSMSMODE_ENUM_NAME : int
{
   MY_LSMSMODE_ENUM
};

#undef X

constexpr const char* lsmsModeToString(LSMSMode e) noexcept
{
#define X(name) case(MY_LSMSMODE_ENUM_NAME::name): return #name;
  switch(e)
  {
    MY_LSMSMODE_ENUM
  }
#undef X
  return "main";
}

#endif
