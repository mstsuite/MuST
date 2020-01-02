#ifndef __IOTK
#define __IOTK

#ifdef LowerCase

#define c_dtsize          c_dtsize
#define c_gopen           c_gopen
#define c_close           c_close
#define c_string_padsize  c_string_padsize
#define c_fseek           c_fseek
#define c_write_double    c_write_double
#define c_write_integer   c_write_integer
#define c_write_string    c_write_string
#define c_read_double     c_read_double
#define c_read_integer    c_read_integer
#define c_read_string     c_read_string

#elif Underscore

#define c_dtsize          c_dtsize_
#define c_gopen           c_gopen_
#define c_close           c_close_
#define c_string_padsize  c_string_padsize_
#define c_fseek           c_fseek_
#define c_write_double    c_write_double_
#define c_write_integer   c_write_integer_
#define c_write_string    c_write_string_
#define c_read_double     c_read_double_
#define c_read_integer    c_read_integer_
#define c_read_string     c_read_string_

#elif DoubleUnderscore

#define c_dtsize          c_dtsize__
#define c_gopen           c_gopen__
#define c_close           c_close__
#define c_string_padsize  c_string_padsize__
#define c_fseek           c_fseek__
#define c_write_double    c_write_double__
#define c_write_integer   c_write_integer__
#define c_write_string    c_write_string__
#define c_read_double     c_read_double__
#define c_read_integer    c_read_integer__
#define c_read_string     c_read_string__

#elif UpperCase

#define c_dtsize          C_DTSIZE
#define c_gopen           C_GOPEN
#define c_close           C_CLOSE
#define c_string_padsize  C_STRING_PADSIZE
#define c_fseek           C_FSEEK
#define c_write_double    C_WRITE_DOUBLE
#define c_write_integer   C_WRITE_INTEGER
#define c_write_string    C_WRITE_STRING
#define c_read_double     C_READ_DOUBLE
#define c_read_integer    C_READ_INTEGER
#define c_read_string     C_READ_STRING

#endif

#endif
