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

#include <wchar.h>

#ifndef RSTRING
#define RSTRING

#ifdef __cplusplus
extern "C"
{
#endif


	struct rui_wstring
	{
		wchar_t *text;	/*<! wchar_t c-string */
		size_t max;	/*<! usable memory allocated to text minus the space for the nul character */
	};

	typedef struct rui_wstring rwstring;

	struct rui_cstring
	{
		char *text;	/*<! char c-string */
		size_t max;	/*<! usable memory allocated to text minus the space for the nul character */
	};

	typedef struct rui_cstring rcstring;


	enum rui_string_error_codes
	{ RS_MEMORY, RS_OK = 1, RS_UNKNOWN };

	typedef enum rui_string_error_codes rstring_code;

/**
Creates a new rwstring by allocating memory and setting up the variables
\param length the allocated size for the new string
\return a pointer to the newly created rwstring
**/
	rwstring *rws_create (size_t length);

/**
Frees the memory allocated to a rwstring
\param rws pointer to a pointer to a rwstring
**/
	void rws_free (rwstring ** rws);

/**
Resizes the maximum memory allocated to the rwstring's wchar_t string
\param rws pointer to the rwstring to be resized
\param length the new maximum memory allocated to rwstring's wchar_t string
\return result
**/
	rstring_code rws_resize (rwstring * rws, size_t length);

/**
Creates a duplicate string of copied
\param copied the string that will be duplicated
\return a duplicate of copied, NULL if the memory allocation failed
**/
	rwstring *rws_duplicate (rwstring * copied);

/** Returns the length of a rwstring
\param rws rwstring to measure
\return the length of the string
**/
	size_t rws_length (rwstring * rws);

/** Copies a rwstring into another rwstring
\param to rwstring where to copy to
\param from rwstring where to copy from
\return result
**/
	rstring_code rws_copyrws (rwstring * to, const rwstring * from);

/** Copies a wchar_t string into a rwstring
\param to rwstring where to copy to
\param from wchar_t string where to copy from
\param length the length to be copied
\return result
**/
	rstring_code rws_copywcs (rwstring * to, const wchar_t * from, const size_t length);

/** Concatenates a rwstring onto the end of a rwstring
\param pre rwstring where to append to
\param pos rtring where to append from
\return result
**/
	rstring_code rws_catrws (rwstring * pre, const rwstring * pos);

/** Concatenates a rcstring onto the end of a rwstring
\param pre rwstring where to append to
\param pos rtring where to append from
\return result
**/
	rstring_code rws_catrcs (rwstring * pre, const rcstring * pos);

/** Concatenates a wchar_t string onto the end of a rwstring
\param pre rwstring where to append to
\param pos wchar_t string where to append from
\param length the length to be copied
\return result
**/
	rstring_code rws_catwcs (rwstring * pre, const wchar_t * pos, const size_t length);

/** Concatenates an UTF-8 c-string onto the end of a rwstring
\param pre rwstring where to append to
\param pos UTF-8 c-string where to append from
\param length the number of characters to be copied
\return result
**/
	rstring_code rws_catcs (rwstring * pre, const char *pos, const size_t length);

/** Concatenates a single wchar_t onto the end of a rwstring
\param pre rwstring where to append to
\param c wchar_t to be appended
\return result
**/
	rstring_code rws_catwc (rwstring * pre, const wchar_t c);

/** Concatenates a single char onto the end of a rwstring
\param pre rwstring where to append to
\param c char to be appended
\return result
**/
	rstring_code rws_catc (rwstring * pre, const char c);

/** Wraps a wchar_t string with a rwstring structure, in order to offer a higher level string handling
\param wcs wchar_t structure which will be wrapped
\return a rwstring structure with the s pointer pointing towards wcs
**/
	rwstring *rws_wrap (wchar_t * wcs);

/** Unwraps a rwstring container from it's wchar_t string. Returns a pointer to the wchar_t string and frees everything else
\param rws the rwstring structure to be unwrapped and freed
\return a wchar_t string containing the text
**/
	wchar_t *rws_unwrap (rwstring * rws);


/**
Creates a new rcstring by allocating memory and setting up the variables
\param length the allocated size for the new string
\return a pointer to the newly created rcstring
**/
	rcstring *rcs_create (size_t length);

/**
Frees the memory allocated to a rcstring
\param rcs pointer to a pointer to a rcstring
**/
	void rcs_free (rcstring ** rcs);

/**
Resizes the maximum memory allocated to the rcstring's char string
\param rcs pointer to the rcstring to be resized
\param length the new maximum memory allocated to rcstring's char string
\return result
**/
	rstring_code rcs_resize (rcstring * rcs, size_t length);

/**
Creates a duplicate string of copied
\param copied the string that will be duplicated
\return a duplicate of copied, NULL if the memory allocation failed
**/
	rcstring *rcs_duplicate (rcstring * copied);

/** Returns the length of a rcstring
\param rcs rcstring to measure
\return the length of the string
**/
	size_t rcs_length (rcstring * rcs);

/** Copies a rcstring into another rcstring
\param to rcstring where to copy to
\param from rcstring where to copy from
\return result
**/
	rstring_code rcs_copyrcs (rcstring * to, const rcstring * from);

/** Copies a char string into a rcstring
\param to rcstring where to copy to
\param from char string where to copy from
\param length the length to be copied
\return result
**/
	rstring_code rcs_copycs (rcstring * to, const char *from, const size_t length);

/** Concatenates a rcstring onto the end of a rcstring
\param pre rcstring where to append to
\param pos rtring where to append from
\return result
**/
	rstring_code rcs_catrcs (rcstring * pre, const rcstring * pos);

/** Concatenates a char string onto the end of a rcstring
\param pre rcstring where to append to
\param pos char string where to append from
\param length the length to be copied
\return result
**/
	rstring_code rcs_catcs (rcstring * pre, const char *pos, const size_t length);

/** Concatenates a wchar_t string onto the end of a rcstring
\param pre rcstring where to append to
\param pos wchar_t string where to append from
\param length the length to be copied
\return result
**/
	rstring_code rcs_catwcs (rcstring * pre, const wchar_t * pos, const size_t length);

/** Concatenates a single char onto the end of a rcstring
\param pre rcstring where to append to
\param c char to be appended
\return result
**/
	rstring_code rcs_catwc (rcstring * pre, const wchar_t c);

/** Concatenates a single char onto the end of a rcstring
\param pre rcstring where to append to
\param c char to be appended
\return result
**/
	rstring_code rcs_catc (rcstring * pre, const char c);

/** Wraps a c-string with a rcstring structure, in order to offer a higher level string handling
\param cs char which will be wrapped
\return a rcstring structure with the s pointer pointing towards cs
**/
	rcstring *rcs_wrap (char *cs);

/** Unwraps a rcstring container from it's char string. Returns a pointer to the char string and frees everything else
\param rcs the rcstring structure to be unwrapped and freed
\return a char string containing the text
**/
	char *rcs_unwrap (rcstring * rcs);


#ifdef __cplusplus
}
#endif
#endif
