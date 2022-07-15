/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2016 M. Oliveira

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

int xc_func_info_get_number(const xc_func_info_type *info)
{
  return info->number;
}

int xc_func_info_get_kind(const xc_func_info_type *info)
{
  return info->kind;
}

char const *xc_func_info_get_name(const xc_func_info_type *info)
{
  return info->name;
}

int xc_func_info_get_family(const xc_func_info_type *info)
{
  return info->family;
}

int xc_func_info_get_flags(const xc_func_info_type *info)
{
  return info->flags;
}

const func_reference_type *xc_func_info_get_references(const xc_func_info_type *info, int number)
{
  assert(number >=0 && number < XC_MAX_REFERENCES);

  if (info->refs[number] == NULL) {
    return NULL;
  } else {
    return info->refs[number];
  }
}

int xc_func_info_get_n_ext_params(const xc_func_info_type *info)
{
  assert(info!=NULL);

  return info->ext_params.n;
}

char const *xc_func_info_get_ext_params_name(const xc_func_info_type *info, int number)
{
  assert(number >=0 && number < info->ext_params.n);

  return info->ext_params.names[number];
}

char const *xc_func_info_get_ext_params_description(const xc_func_info_type *info, int number)
{
  assert(number >=0 && number < info->ext_params.n);

  return info->ext_params.descriptions[number];
}

double xc_func_info_get_ext_params_default_value(const xc_func_info_type *info, int number)
{
  assert(number >=0 && number < info->ext_params.n);

  return info->ext_params.values[number];
}
