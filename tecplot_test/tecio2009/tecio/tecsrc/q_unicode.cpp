/*
 * NOTICE and LICENSE for Tecplot Input/Output Library (TecIO) - OpenFOAM
 *
 * Copyright (C) 1988-2009 Tecplot, Inc.  All rights reserved worldwide.
 *
 * Tecplot hereby grants OpenCFD limited authority to distribute without
 * alteration the source code to the Tecplot Input/Output library, known 
 * as TecIO, as part of its distribution of OpenFOAM and the 
 * OpenFOAM_to_Tecplot converter.  Users of this converter are also hereby
 * granted access to the TecIO source code, and may redistribute it for the
 * purpose of maintaining the converter.  However, no authority is granted
 * to alter the TecIO source code in any form or manner.
 *
 * This limited grant of distribution does not supersede Tecplot, Inc.'s 
 * copyright in TecIO.  Contact Tecplot, Inc. for further information.
 * 
 * Tecplot, Inc.
 * 3535 Factoria Blvd, Ste. 550
 * Bellevue, WA 98006, USA
 * Phone: +1 425 653 1200
 * http://www.tecplot.com/
 *
 */
#include "stdafx.h"
#include "MASTER.h"

#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2008 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#define Q_UNICODEMODULE

#include "GLOBAL.h"
#include "TASSERT.h"

#if !defined TECPLOTKERNEL
#include "TranslatedString.h"
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#include "ALLOC.h"

#include "Q_UNICODE.h"

using namespace std;

namespace tecplot
{
namespace strutil
{

typedef std::map<std::string, char *>      EnvStringPoolMap_t;
static EnvStringPoolMap_t       mapEnvStringPool;


#if defined MSWIN


string WStringToString(wstring str)
{
    REQUIRE("str is any wide string");
    string Result = WideCharToUtf8(str.c_str());

    ENSURE("Result is any string");
    return Result;
}

wstring StringToWString(string str)
{
    REQUIRE("str is any string");

    wstring Result = Utf8ToWideChar(str.c_str());

    ENSURE("Result is any string");
    return Result;
}
#endif

/************************************************
 * Utf8Api
 ************************************************/
#define VALID_CODE_PAGE(cp) \
  ( (cp) == 932 || (cp) == CP_UTF8 || (cp) == CP_ACP || (cp) == CP_OEMCP || (cp) == CP_THREAD_ACP )


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined MSWIN && !defined ENGINE
#endif
#if defined MSWIN
#endif
#endif


#if defined MSWIN

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

std::string  WideCharToMultiByte(const wchar_t* WideString,
                                 unsigned int       CodePage)
{
    REQUIRE(VALID_REF(WideString));
    REQUIRE(VALID_CODE_PAGE(CodePage));

    string  strResult;
    wstring wString(WideString);


    if (wString.length() > 0)
    {
        size_t nLen =
            (size_t) ::WideCharToMultiByte(CodePage,
                                           0,
                                           wString.c_str(),
                                           -1,
                                           NULL,
                                           0,
                                           NULL,
                                           NULL);
        if (nLen > 0)
        {
            char *pBuffer = ALLOC_ARRAY(nLen, char, "pBuffer");

            VERIFY(::WideCharToMultiByte(CodePage,
                                         0,
                                         WideString,
                                         (int)(wString.length() + 1),
                                         pBuffer,
                                         (int)nLen,
                                         NULL,
                                         NULL) != 0);

            strResult = pBuffer;
            FREE_ARRAY(pBuffer, "pBuffer");

        }
        else
        {
            // this should never be an error
            CHECK(FALSE);
        }
    }
    else
    {
        // output 'str' remains empty
    }

    ENSURE("strResult is a valid STL string");
    return strResult;


}

wstring MultiByteToWideChar(const char     *UTF8String,
                            unsigned int    CodePage)
{
    REQUIRE(VALID_REF(UTF8String));
    REQUIRE(VALID_CODE_PAGE(CodePage));

    wstring strResult;
    string  UTF8str(UTF8String);

    size_t wLen;

    if (UTF8str.length() > 0)
    {
        wLen =
            (size_t) ::MultiByteToWideChar(CodePage,
                                           0,
                                           UTF8str.c_str(),
                                           -1,
                                           NULL,
                                           0);
        if (wLen > 0)
        {
            wchar_t *wBuffer = ALLOC_ARRAY(wLen + 1, wchar_t, "wBuffer");
            VERIFY(::MultiByteToWideChar(CodePage,
                                         0,
                                         UTF8str.c_str(),
                                         (int)(UTF8str.length() + 1),
                                         wBuffer,
                                         (int)wLen) != 0);

            strResult = wBuffer;
            FREE_ARRAY(wBuffer, "wBuffer");

        }
        else
        {
            CHECK(FALSE); // We should never get an error here
        }
    }
    else
    {
        // strResult is left empty
    }

    ENSURE("strResult is a valid CString");

    wstring strRet(strResult);
    return strRet;

}
#endif



#if defined MSWIN
std::string WideCharToUtf8(const wchar_t *str)
{
    REQUIRE(VALID_REF(str)); /* really cannot be NULL - 2007-10-22 CAM/DTO */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    UINT CodePage = CP_ACP;

    string Result = "";

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    Result = WideCharToMultiByte(str, CodePage);

    ENSURE("Result is any string");
    return Result;
}

wstring Utf8ToWideChar(const char *str)
{
    REQUIRE(VALID_REF(str)); /* really cannot be NULL - 2007-10-22 CAM/DTO */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    UINT CodePage = CP_ACP;
    wstring Result;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    Result = MultiByteToWideChar(str, CodePage);

    ENSURE("Result is any string");
    return Result;
}
#endif


Boolean_t IsValidUtf8LeadByte(Byte_t uch)
{
    REQUIRE("uch is any byte");
    Boolean_t Result =
        uch <= 0x7F                   ||
        (uch >= 0xC0 && uch <= 0xDF) ||
        (uch >= 0xE0 && uch <= 0xEF) ||
        (uch >= 0xF0 && uch <= 0xF4);

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}

Boolean_t IsValidUtf8ContinuingByte(Byte_t uch)
{
    REQUIRE("uch is any char");

    Boolean_t Result =
        (uch >= 0x80 && uch <= 0xBF);

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}

Boolean_t IsValidUtf8Byte(Byte_t uch)
{
    REQUIRE("uch is any char");
    Boolean_t Result =
        IsValidUtf8LeadByte(uch)        ||
        IsValidUtf8ContinuingByte(uch);

    REQUIRE(VALID_BOOLEAN(Result));
    return Result;
}


Boolean_t ShouldConvertWideStringToUtf8String(const wchar_t *str)
{
    Boolean_t Result = FALSE;

#if defined MSWIN && defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* MSWIN and TECPLOTKERNEL */

    ENSURE(VALID_BOOLEAN(Result));
    return Result;

}

Boolean_t IsValidUtf8String(const char *str)
{
    Boolean_t IsValid = TRUE;
    REQUIRE(VALID_REF(str));

#if defined MSWIN
    size_t len                      = strlen(str);
    Boolean_t InUtf8Sequence        = FALSE;
    int       Utf8SequenceCount     = 0;
    int       MaxUtf8SequenceCount  = 0;

    /* we want to process the final \0, so go to <= len */

    for (size_t ii = 0; IsValid && ii <= len; ii++)
    {
        Byte_t uch = (Byte_t)str[ii];

        if (uch <= 0x7F)
        {
            /* This must be the end of a sequence,
               so the sequence count must match
               the max sequence count */

            InUtf8Sequence        = FALSE;
            IsValid               = (Utf8SequenceCount == MaxUtf8SequenceCount);
            Utf8SequenceCount     = 0;
            MaxUtf8SequenceCount  = 0;
        }
        else if (uch >= 0x80 && uch <= 0xBF)
        {
            /* Continuing byte in a multi byte sequence */
            if (InUtf8Sequence)
            {
                Utf8SequenceCount++;
            }
            else
            {
                IsValid = FALSE;
            }

        }
        else if (uch >= 0xC0 && uch <= 0xDF)
        {
            /* Lead byte of 000080-0007FF */
            IsValid               = (Utf8SequenceCount == MaxUtf8SequenceCount);
            InUtf8Sequence        = TRUE;
            Utf8SequenceCount     = 0;
            MaxUtf8SequenceCount  = 1;
        }
        else if (uch >= 0xE0 && uch <= 0xEF)
        {
            /* Lead byte of 000800-00FFFF */
            IsValid               = (Utf8SequenceCount == MaxUtf8SequenceCount);
            InUtf8Sequence        = TRUE;
            Utf8SequenceCount     = 0;
            MaxUtf8SequenceCount  = 2;
        }
        else if (uch >= 0xF0 && uch <= 0xF4)
        {
            /* Lead byte of 010000-10FFFF */
            IsValid               = (Utf8SequenceCount == MaxUtf8SequenceCount);
            Utf8SequenceCount     = 0;
            InUtf8Sequence        = TRUE;
            MaxUtf8SequenceCount  = 3;
        }

        else
        {
            /* Invalid Utf 8 */
            IsValid = FALSE;
        }
    }
#endif

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}


/**
 */
Boolean_t IsNullOrZeroLengthString(const char *str)
{
    REQUIRE(VALID_REF_OR_NULL(str));

    Boolean_t Result = (str == NULL || strlen(str) == 0);

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}

/**
 */
Boolean_t IsNullOrZeroLengthString(TranslatedString TS)
{
    REQUIRE(TS.isValid());
    return TS.isNullOrZeroLength();
}


}
}

#if defined MSWIN && TECPLOTKERNEL && (!defined NO_ASSERTS || defined CHECKED_BUILD)
/* Keeping Trace out of the release builds
   will verify for us that it has been optimized away.

   See the definition of TRACE in MASTER.h for
   more information... */
void MSWinTrace(const char *Format, ...)
{
    REQUIRE(VALID_REF(Format));

    const int BufferSize = 512; /* Only print the first 512 characers */
    va_list Arguments;

    /* Don't use ALLOC_ARRAY here */
    char *buffer = new char[BufferSize];
    memset(buffer, 0, BufferSize);

    va_start(Arguments, Format);
    _vsnprintf(buffer, BufferSize - 1, Format, Arguments);
    va_end(Arguments);

    ::OutputDebugStringA(buffer);

    delete [] buffer;
}



#endif
