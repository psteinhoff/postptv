#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <windows.h>
#include "typedefs.h"
#include "globals.h"

#include <tools.h>

#pragma once

#include <tchar.h>
#include <conio.h>

int PostToParaview_c( int argc, const char** argv, char base_path[256],int flag4export[4]);