/*
 *      DiversUtils - Functions to measure diversity
 *
 * Copyright (c) 2025  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdint.h>
#include <string.h>

#include "case.h"

void lcase_ascii(char * restrict const s){
    int64_t i = 0;
    while(s[i]){
        if('A' <= s[i] && s[i] <= 'Z'){s[i] += ('a' - 'A');}
        i++;
    }
}

void lcase_utf8_french(void * restrict const s){
    int64_t i = 0;
    /*
    while(s[i]){
        if('A' <= s[i] && s[i] <= 'Z'){s[i] += ('a' - 'A'); i++;}
        else if(s[i] <= 0x7f){i++;}
        else if(strcmp(&s[i], "À") == 0){strcpy(&s[i], "à"); i += strlen("à");}
        else if(strcmp(&s[i], "Â") == 0){strcpy(&s[i], "â"); i += strlen("â");}
        else if(strcmp(&s[i], "Ä") == 0){strcpy(&s[i], "ä"); i += strlen("ä");}
        else if(strcmp(&s[i], "É") == 0){strcpy(&s[i], "é"); i += strlen("é");}
        else if(strcmp(&s[i], "È") == 0){strcpy(&s[i], "è"); i += strlen("è");}
        else if(strcmp(&s[i], "Ê") == 0){strcpy(&s[i], "ê"); i += strlen("ê");}
        else if(strcmp(&s[i], "Ë") == 0){strcpy(&s[i], "ë"); i += strlen("ë");}
        else if(strcmp(&s[i], "Ì") == 0){strcpy(&s[i], "ì"); i += strlen("ì");}
        else if(strcmp(&s[i], "Î") == 0){strcpy(&s[i], "î"); i += strlen("î");}
        else if(strcmp(&s[i], "Ï") == 0){strcpy(&s[i], "ï"); i += strlen("ï");}
        else if(strcmp(&s[i], "Ò") == 0){strcpy(&s[i], "ò"); i += strlen("ò");}
        else if(strcmp(&s[i], "Ô") == 0){strcpy(&s[i], "ô"); i += strlen("ô");}
        else if(strcmp(&s[i], "Ö") == 0){strcpy(&s[i], "ö"); i += strlen("ö");}
        else if(strcmp(&s[i], "Ù") == 0){strcpy(&s[i], "ù"); i += strlen("ù");}
        else if(strcmp(&s[i], "Û") == 0){strcpy(&s[i], "û"); i += strlen("û");}
        else if(strcmp(&s[i], "Ü") == 0){strcpy(&s[i], "ü"); i += strlen("ü");}
        else {i++;}
    }
    */
    while(((unsigned char *) s)[i]){
        if('A' <= ((unsigned char *) s)[i] && ((unsigned char *) s)[i] <= 'Z'){((unsigned char *) s)[i] += ('a' - 'A'); i++;}
        else if(((unsigned char *) s)[i] <= 0x7f){i++;}
        else if(strcmp(&((char *)s)[i], "Ç") == 0){strcpy(&((char *)s)[i], "ç"); i += strlen("ç");}
        else if(strcmp(&((char *)s)[i], "À") == 0){strcpy(&((char *)s)[i], "à"); i += strlen("à");}
        else if(strcmp(&((char *)s)[i], "Â") == 0){strcpy(&((char *)s)[i], "â"); i += strlen("â");}
        else if(strcmp(&((char *)s)[i], "Ä") == 0){strcpy(&((char *)s)[i], "ä"); i += strlen("ä");}
        else if(strcmp(&((char *)s)[i], "É") == 0){strcpy(&((char *)s)[i], "é"); i += strlen("é");}
        else if(strcmp(&((char *)s)[i], "È") == 0){strcpy(&((char *)s)[i], "è"); i += strlen("è");}
        else if(strcmp(&((char *)s)[i], "Ê") == 0){strcpy(&((char *)s)[i], "ê"); i += strlen("ê");}
        else if(strcmp(&((char *)s)[i], "Ë") == 0){strcpy(&((char *)s)[i], "ë"); i += strlen("ë");}
        else if(strcmp(&((char *)s)[i], "Ì") == 0){strcpy(&((char *)s)[i], "ì"); i += strlen("ì");}
        else if(strcmp(&((char *)s)[i], "Î") == 0){strcpy(&((char *)s)[i], "î"); i += strlen("î");}
        else if(strcmp(&((char *)s)[i], "Ï") == 0){strcpy(&((char *)s)[i], "ï"); i += strlen("ï");}
        else if(strcmp(&((char *)s)[i], "Ò") == 0){strcpy(&((char *)s)[i], "ò"); i += strlen("ò");}
        else if(strcmp(&((char *)s)[i], "Ô") == 0){strcpy(&((char *)s)[i], "ô"); i += strlen("ô");}
        else if(strcmp(&((char *)s)[i], "Ö") == 0){strcpy(&((char *)s)[i], "ö"); i += strlen("ö");}
        else if(strcmp(&((char *)s)[i], "Ù") == 0){strcpy(&((char *)s)[i], "ù"); i += strlen("ù");}
        else if(strcmp(&((char *)s)[i], "Û") == 0){strcpy(&((char *)s)[i], "û"); i += strlen("û");}
        else if(strcmp(&((char *)s)[i], "Ü") == 0){strcpy(&((char *)s)[i], "ü"); i += strlen("ü");}
        else {i++;}
    }
}
