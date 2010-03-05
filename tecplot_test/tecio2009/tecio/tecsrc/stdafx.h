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
#if !defined STDAFX_H_
# define STDAFX_H_

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2008 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
/*
 * stdafx.h : include file for standard system include files,
 *  or project specific include files that are used frequently, but
 *      are changed infrequently
 * Used for Windows only
 */
#if defined _WIN32

/*
 * Set NDEBUG before including "custafx.h" since that file may
 * use NDEBUG.  (In fact, for SmartHeap builds this is the case.)
 * CAM 2007-04-11
 *
 * Previous comment: "Note that _DEBUG is defined by the Windows compiler
 * if any of the multi-threaded DLL runtime libraries are used."
 */
#if !defined _DEBUG
#if !defined NDEBUG
#define NDEBUG
#endif
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#if !defined MSWIN
#define MSWIN
#endif

#define ENGLISH_ONLY // remove to support non-english dll's

#if !defined WINVER
#define WINVER 0x0500
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined CHECKED_BUILD || defined _DEBUG && !defined COREAPI
#if defined _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
#endif
#else
#endif
#endif /* TECPLOTKERNEL */

/*  Windows builds are UNICODE */
#pragma warning(disable : 4786) /* truncated identifiers in debug symbol table. */
#pragma warning(disable : 4996) /* deprecated functions */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined UNICODE
#endif
#if !defined _AFXDLL
#endif
#if defined _M_IX86
#else
#endif
#if defined _WIN64
#if !defined _M_IA64 && !defined _M_AMD64
#endif
#endif
#if !defined MSWIN
#endif
#if !defined THREED
#endif
#ifndef _AFX_NO_AFXCMN_SUPPORT
#endif /* _AFX_NO_AFXCMN_SUPPORT */
#ifndef _AFX
#endif
#else /* !TECPLOTKERNEL */
#define AfxIsValidAddress(ptr,bb) ((ptr)!=NULL)
#endif

/* 64-bit adjustments */
#if defined _M_IA64 || defined _M_AMD64
#define WININT  INT_PTR
#define WINUINT UINT_PTR
#else
#define WININT  int
#define WINUINT UINT
#endif

#define WINCALLBACK CALLBACK

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined (NDEBUG)
#else
#endif
#endif /* TECPLOTKERNEL */
#endif /* _WIN32 */


#endif /* STDAFX_H_ */
