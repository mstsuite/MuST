/***************************************************************************
 *   Copyright (C) 2007 by Rui Maciel   *
 *   rui.maciel@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "rstring.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <wchar.h>


#define RSTRING_INCSTEP 3


rwstring *
rws_create (size_t length)
{
	rwstring *rws;
	rws = malloc (sizeof (rwstring));	/* allocates memory for a struct rwstring */
	if (rws == NULL)
		return NULL;

	rws->text = calloc (length + 1, sizeof (wchar_t));	/* text[max] reserved for L'\0' */
	if (rws->text == NULL)
	{
		free (rws);
		return NULL;
	}

	rws->max = length;
	return rws;
}


void
rws_free (rwstring ** rws)
{
	assert (rws != NULL);
	if (*rws != NULL)
	{
		if ((*rws)->text != NULL)
		{
			free ((*rws)->text);
			(*rws)->text = NULL;
		}
		free (*rws);
		*rws = NULL;
	}

}


rstring_code
rws_resize (rwstring * rws, size_t length)
{
	wchar_t *temp;
	assert (rws != NULL);

	temp = realloc (rws->text, sizeof (wchar_t) * (length + 1));	/* length plus L'\0' */
	if (temp == NULL)
	{
		free (rws);
		return RS_MEMORY;
	}
	rws->text = temp;
	rws->max = length;
	rws->text[rws->max] = L'\0';
	return RS_OK;
}


rwstring *
rws_duplicate (rwstring * copied)
{
	rwstring *copy;
	assert (copied != NULL);

	copy = rws_create (copy->max);
	if (copy == NULL)
		return NULL;

	wcsncpy (copy->text, copied->text, wcslen (copied->text));
	return copy;
}


size_t
rws_length (rwstring * rws)
{
	assert (rws != NULL);
	return wcslen (rws->text);
}


rstring_code
rws_copyrws (rwstring * to, const rwstring * from)
{
	size_t from_length;
	assert (from != NULL);
	assert (to != NULL);

	from_length = wcslen (from->text);
	/*TODO implement intelligent memory allocation */
	if (to->max < from_length)
	{
		if (rws_resize (to, from_length + 1) != RS_OK)
			return RS_MEMORY;
	}
	wcsncpy (to->text, from->text, from_length);
	to->text[from_length] = L'\0';
	return RS_OK;
}


rstring_code
rws_copywcs (rwstring * to, const wchar_t * from, const size_t length)
{
	assert (from != NULL);
	assert (to != NULL);

	/*TODO implement intelligent memory allocation */
	if (to->max < length)
	{
		if (rws_resize (to, length + 1) != RS_OK)
			return RS_MEMORY;
	}
	wcsncpy (to->text, from, length);
	to->text[length] = L'\0';
	return RS_OK;
}


rstring_code
rws_catrws (rwstring * pre, const rwstring * pos)
{
	size_t pre_length, pos_length;
	assert (pre != NULL);
	assert (pos != NULL);

	pre_length = wcslen (pre->text);
	pos_length = wcslen (pos->text);
	if (pos == NULL)
		return RS_OK;

	if (pre->max < pre_length + pos_length + 1)
	{
		if (rws_resize (pre, pre_length + pos_length + 1) != RS_OK)
			return RS_MEMORY;
	}
	wcsncpy (pre->text + pre_length, pos->text, pos_length);
	pre->text[pre_length + pos_length] = L'\0';
	return RS_OK;
}


rstring_code
rws_catrcs (rwstring * pre, const rcstring * pos)
{
	size_t utf8pos;
	wchar_t wc;
	char i;			/* static loop counter */

	assert (pre != NULL);
	assert (pos != NULL);

	/* starting the conversion */
	utf8pos = 0;

	while (utf8pos < strlen (pos->text))
	{
		if ((pos->text[utf8pos] & 0x80) == 0)
		{
			if (rws_catc (pre, pos->text[utf8pos++]) != RS_OK)
			{
				return RS_MEMORY;
			}
		}
		else if ((pos->text[utf8pos] & 0xE0) == 0xC0)
		{
			wc = (pos->text[utf8pos++] & 0x1F) << 6;
			if ((pos->text[utf8pos] & 0xC0) == 0x80)
			{
				wc |= (pos->text[utf8pos++] & 63);

			}
			else
			{	/* Invalid utf8 string */
				return RS_UNKNOWN;	/* malformed UTF8 string */
			}
			if (rws_catwc (pre, wc) != RS_OK)
			{
				return RS_MEMORY;
			}
		}
		else if ((pos->text[utf8pos] & 0xF0) == 0xE0)
		{
			wc = pos->text[utf8pos++] & 0xF << 12;
			for (i = 1; i >= 0; i--)
			{
				if ((pos->text[utf8pos] & 0xC0) == 0x80)
				{
					wc &= pos->text[utf8pos++] & 0x3F << i;
				}
				else
				{	/* Invalid utf8 string */
					return RS_UNKNOWN;
				}
			}
			if (rws_catwc (pre, wc) != RS_OK)
			{
				return RS_MEMORY;
			}

		}
		else if ((pos->text[utf8pos] & 0xF8) == 0xF0)
		{
			wc = pos->text[utf8pos++] & 0xF << 12;
			for (i = 2; i >= 0; i--)
			{
				if ((pos->text[utf8pos] & 0xC0) == 0x80)
				{
					wc &= pos->text[utf8pos++] & 0x3F << i;
				}
				else
				{	/* Invalid utf8 string */
					return RS_UNKNOWN;

				}
			}
			if (rws_catwc (pre, wc) != RS_OK)
			{
				return RS_MEMORY;
			}

		}
		else if ((pos->text[utf8pos] & 0xFC) == 0xF8)
		{
			for (i = 3; i >= 0; i--)
			{
				if ((pos->text[utf8pos] & 0xC0) == 0x80)
				{
					wc &= pos->text[utf8pos++] & 0x3F << i;
				}
				else
				{	/* Invalid utf8 string */
					return RS_UNKNOWN;

				}
			}
			if (rws_catwc (pre, wc) != RS_OK)
			{
				return RS_MEMORY;
			}
		}
		else if ((pos->text[utf8pos] & 0xFE) == 0xFC)
		{
			for (i = 4; i >= 0; i--)
			{
				if ((pos->text[utf8pos] & 0xC0) == 0x80)
				{
					wc &= pos->text[utf8pos++] & 0x3F << i;
				}
				else
				{	/* Invalid utf8 string */
					return RS_UNKNOWN;
				}
			}
			if (rws_catwc (pre, wc) != RS_OK)
			{
				return RS_MEMORY;
			}
		}
		else		/* Invalid utf8 string */
		{
			return RS_UNKNOWN;
		}
	}

	return RS_OK;
}


rstring_code
rws_catwcs (rwstring * pre, const wchar_t * pos, const size_t length)
{
	size_t pre_length;

	assert (pre != NULL);
	assert (pos != NULL);

	pre_length = wcslen (pre->text);

	if (pre->max < pre_length + length)
	{
		if (rws_resize (pre, pre_length + length + 1) != RS_OK)
			return RS_MEMORY;
	}
	wcsncpy (pre->text + pre_length, pos, length);
	pre->text[pre_length + length] = L'\0';
	return RS_OK;
}


rstring_code
rws_catcs (rwstring * pre, const char *pos, const size_t length)
{
	size_t pre_length;

	assert (pre != NULL);
	assert (pos != NULL);

	pre_length = wcslen (pre->text);

	if (pre->max < pre_length + length)
	{
		if (rws_resize (pre, pre_length + length + 1) != RS_OK)
			return RS_MEMORY;
	}
	mbsrtowcs (pre->text + pre_length, &pos, length, NULL);
	pre->text[pre_length + length] = L'\0';
	return RS_OK;
}


rstring_code
rws_catwc (rwstring * pre, const wchar_t c)
{
	size_t pre_length;

	assert (pre != NULL);

	pre_length = wcslen (pre->text);
	if (pre->max <= pre_length)
	{
		pre->max += RSTRING_INCSTEP;
		if (rws_resize (pre, pre->max) != RS_OK)
			return RS_MEMORY;
	}
	pre->text[pre_length] = c;
	pre->text[pre_length + 1] = L'\0';
	return RS_OK;
}


rstring_code
rws_catc (rwstring * pre, const char c)
{
	wchar_t newc;

	assert (pre != NULL);
	mbtowc (&newc, &c, 1);
	return rws_catwc (pre, newc);
}


rwstring *
rws_wrap (wchar_t * wcs)
{
	rwstring *wrapper;

	assert (wcs != NULL);
	wrapper = malloc (sizeof (rwstring));
	if (wrapper == NULL)
		return NULL;
	wrapper->max = wcslen (wcs);
	wrapper->text = wcs;
	return wrapper;
}


wchar_t *
rws_unwrap (rwstring * rws)
{
	wchar_t *out;
	assert (rws != NULL);
	if (rws->text == NULL)
		out = NULL;
	else
		out = realloc (rws->text, sizeof (wchar_t) * (wcslen (rws->text) + 1));
	free (rws);
	return out;
}


rcstring *
rcs_create (size_t length)
{
	rcstring *rcs;
	rcs = malloc (sizeof (rcstring));	/* allocates memory for a struct rcstring */
	if (rcs == NULL)
		return NULL;

	rcs->max = length;

	rcs->text = calloc (rcs->max + 1, sizeof (char));
	if (rcs->text == NULL)
	{
		free (rcs);
		return NULL;
	}

	return rcs;
}


void
rcs_free (rcstring ** rcs)
{
	assert (rcs != NULL);
	if (*rcs != NULL)
	{
		if ((*rcs)->text != NULL)
		{
			free ((*rcs)->text);
			(*rcs)->text = NULL;
		}
		free (*rcs);
		*rcs = NULL;
	}

}


rstring_code
rcs_resize (rcstring * rcs, size_t length)
{
	char *temp;
	assert (rcs != NULL);

	temp = realloc (rcs->text, sizeof (char) * (length + 1));	/* length plus L'\0' */
	if (temp == NULL)
	{
		free (rcs);
		return RS_MEMORY;
	}
	rcs->text = temp;
	rcs->max = length;
	rcs->text[rcs->max] = L'\0';
	return RS_OK;
}


rcstring *
rcs_duplicate (rcstring * copied)
{
	rcstring *copy;
	assert (copied != NULL);

	copy = rcs_create (copy->max);
	if (copy == NULL)
		return NULL;

	strncpy (copy->text, copied->text, strlen (copied->text));
	return copy;
}


size_t
rcs_length (rcstring * rcs)
{
	/*TODO account for UTF8 */
	assert (rcs != NULL);
	return strlen (rcs->text);
}


rstring_code
rcs_copyrcs (rcstring * to, const rcstring * from)
{
	size_t from_length;
	assert (from != NULL);
	assert (to != NULL);

	from_length = strlen (from->text);
	/*TODO implement intelligent memory allocation */
	if (to->max < from_length)
	{
		if (rcs_resize (to, from_length + 1) != RS_OK)
			return RS_MEMORY;
	}
	strncpy (to->text, from->text, from_length);
	to->text[from_length] = '\0';
	return RS_OK;
}


rstring_code
rcs_copycs (rcstring * to, const char *from, const size_t length)
{
	assert (to != NULL);

	if (from == NULL)
		return RS_OK;

	/*TODO implement intelligent memory allocation */
	if (to->max < length)
	{
		if (rcs_resize (to, length + 1) != RS_OK)
			return RS_MEMORY;
	}
	strncpy (to->text, from, length);
	to->text[length] = '\0';
	return RS_OK;
}


rstring_code
rcs_catrcs (rcstring * pre, const rcstring * pos)
{
	size_t pre_length, pos_length;
	assert (pre != NULL);
	assert (pos != NULL);

	pre_length = strlen (pre->text);
	pos_length = strlen (pos->text);

	if (pre->max < pre_length + pos_length + 1)
	{
		if (rcs_resize (pre, pre_length + pos_length + 1) != RS_OK)
			return RS_MEMORY;
	}
	strncpy (pre->text + pre_length, pos->text, pos_length);
	pre->text[pre_length + pos_length] = '\0';
	return RS_OK;
}


rstring_code
rcs_catcs (rcstring * pre, const char *pos, const size_t length)
{
	size_t pre_length;

	assert (pre != NULL);
	assert (pos != NULL);

	pre_length = strlen (pre->text);

	if (pre->max < pre_length + length)
	{
		if (rcs_resize (pre, pre_length + length + 5) != RS_OK)
			return RS_MEMORY;
	}
	strncpy (pre->text + pre_length, pos, length);
	pre->text[pre_length + length] = '\0';
	return RS_OK;
}


rstring_code
rcs_catwcs (rcstring * pre, const wchar_t * pos, const size_t length)
{
	size_t i;
	rstring_code code;
	assert (pre != NULL);
	assert (pos != NULL);

	for (i = 0; i < length; i++)
	{
		code = rcs_catwc (pre, pos[i]);
		if (code != RS_OK)
		{
			return code;
		}
	}

	return RS_OK;
}


rstring_code
rcs_catwc (rcstring * pre, const wchar_t wc)
{
	assert (pre != NULL);
	/* convert wc to multi-byte UTF8 string and append the product */

	if (wc <= 0x7F)
	{
		rcs_catc (pre, wc);
	}
	else if (wc <= 0x7FF)
	{
		rcs_catc (pre, (wc >> 6) | 192);
		rcs_catc (pre, (wc & 63) | 128);
	}
	else if (wc <= 0xFFFF)
	{
		rcs_catc (pre, wc >> 12 | 224);
		rcs_catc (pre, (wc >> 6 & 63) | 128);
		rcs_catc (pre, (wc & 63) | 128);
	}
	else if (wc <= 0x1FFFFF)
	{
		rcs_catc (pre, wc >> 18 | 240);
		rcs_catc (pre, (wc >> 12 & 63) | 128);
		rcs_catc (pre, (wc >> 6 & 63) | 128);
		rcs_catc (pre, (wc & 63) | 128);
	}
	else if (wc <= 0x3FFFFFF)
	{
		rcs_catc (pre, wc >> 24 | 248);
		rcs_catc (pre, (wc >> 18 & 63) | 128);
		rcs_catc (pre, (wc >> 12 & 63) | 128);
		rcs_catc (pre, (wc >> 6 & 63) | 128);
		rcs_catc (pre, (wc & 63) | 128);
	}
	else if (wc <= 0x7FFFFFFF)
	{
		rcs_catc (pre, wc >> 30 | 252);
		rcs_catc (pre, (wc >> 24 & 63) | 128);
		rcs_catc (pre, (wc >> 18 & 63) | 128);
		rcs_catc (pre, (wc >> 12 & 63) | 128);
		rcs_catc (pre, (wc >> 6 & 63) | 128);
		rcs_catc (pre, (wc & 63) | 128);
	}
	return RS_OK;
}


rstring_code
rcs_catc (rcstring * pre, const char c)
{
	size_t pre_length;

	assert (pre != NULL);

	pre_length = strlen (pre->text);
	if (pre->max <= pre_length)
	{
		pre->max += RSTRING_INCSTEP;
		if (rcs_resize (pre, pre->max) != RS_OK)
			return RS_MEMORY;
	}
	pre->text[pre_length] = c;
	pre->text[pre_length + 1] = '\0';
	return RS_OK;
}


rcstring *
rcs_wrap (char *cs)
{
	rcstring *wrapper;
	assert (cs != NULL);

	wrapper = malloc (sizeof (rcstring));
	if (wrapper == NULL)
		return NULL;
	wrapper->max = strlen (cs);
	wrapper->text = cs;
	return wrapper;
}


char *
rcs_unwrap (rcstring * rcs)
{
	char *out;
	assert (rcs != NULL);

	if (rcs->text == NULL)
		out = NULL;
	else
		out = realloc (rcs->text, sizeof (char) * (strlen (rcs->text) + 1));

	free (rcs);
	return out;
}
