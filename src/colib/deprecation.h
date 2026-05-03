// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef DEPRECATED
#if __GNUC__>3
#define DEPRECATED __attribute__((deprecated))
#define PRIVATE __attribute__((deprecated))
#else
#define DEPRECATED
#define PRIVATE
#endif
#endif
