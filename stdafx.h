// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#else
#include <sys/time.h>
#endif
// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <stdio.h>
// TODO: reference additional headers your program requires here
#include <GL/glut.h>

#pragma warning(disable: 4305) /*truncation from 'double' to 'float'*/
#pragma warning(disable: 4244) /*conversion from 'xxx' to 'yyy', possible loss of data*/
