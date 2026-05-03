// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef debugf_h_
#define debugf_h_

namespace {
    void iprintf(FILE *stream,int depth,const char *fmt,...) {
        fprintf(stream,"%*s",depth,"");
        va_list args;
        va_start(args,fmt);
        vfprintf(stream,fmt,args);
        va_end(args);
    }

    bool strflag(const char *list,const char *key) {
        while(*list) {
            // at the beginning of a potential match here
            bool match = true;
            const char *p = key;
            while(*list && *list!=',' && *p) {
                if(*list++!=*p++) {
                    match = false;
                    break;
                }
            }
            // check that the match ended with a separator
            if(match&&(!*list||*list==',')&&!*p)
                return true;
            // advance to the next separator
            while(*list && *list!=',')
                list++;
            if(*list==',')
                list++;
        }
        return false;
    }

    bool debug(const char *which) {
        const char *env_always = getenv("debug_always");
        if(!env_always) env_always = "info,warn,error,fixme";
        const char *env = getenv("debug");
        if(!env) env = "";
        return (strflag(env_always,which) || strflag(env,which));
    }

    void debugf(const char *which,const char *fmt,...) {
        if(!debug(which)) return;
        va_list args;
        fprintf(stderr,"[%s] ",which);
        va_start(args,fmt);
        vfprintf(stderr,fmt,args);
        va_end(args);
    }
}

#endif
