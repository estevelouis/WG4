/*
 *      DiversUtils - Functions to measure diversity
 *
 * Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
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

#define DEFINE_UNICODE_CONSTANTS
#include "unicode/utf8.h"
#include "unicode/unicode.h"

int8_t utf8_get_length(const char * const s){
    if((s[0] & 0xf8) == 0xf0){return 4;}      // if((s[0] & 0b11111000) == 0b11110000){return 4;} // starts with 11110
    else if((s[0] & 0xf0) == 0xe0){return 3;} // else if((s[0] & 0b11110000) == 0b11100000){return 3;} // starts with 1110
    else if((s[0] & 0xe0) == 0xc0){return 2;} // else if((s[0] & 0b11100000) == 0b11000000){return 2;} // starts with 110
    else if((s[0] & 0x80) == 0x00){return 1;} // else if((s[0] & 0b10000000) == 0b00000000){return 1;} // starts with 0
    else {return -1;} // not well-formed
}

int32_t utf8_is_well_formed(const char * const s, int8_t * const len_p){
    const int8_t len = utf8_get_length(s);

    *len_p = len;

    if(len == -1){return 0;} // not well-formed

    int8_t i;
    for(i = 1 ; i < len ; i++){
        if((s[i] & 0xc0) != 0x80){return 0;} // if((s[i] & 0b11000000) != 0b10000000){return 0;} // not well-formed
    }

    return 1; // well-formed
}

unicode_t utf8_to_unicode(const char * const s, int8_t * const len_p){
    unicode_t unicode = 0;
    int8_t len = utf8_get_length(s);

    *len_p = len;

    if(len == -1){return -1;}

    if(len == 1){
        unicode = (unicode_t) ((uint8_t) s[0]);
    } else if(len == 2){
        unicode = ((unicode_t) ((uint8_t) s[0] & 0x1f)) << 6; // unicode = ((unicode_t) ((uint8_t) s[0] & 0b00011111)) << 6;
        unicode += ((unicode_t) ((uint8_t) s[1] & 0x3f)); // unicode += ((unicode_t) ((uint8_t) s[1] & 0b00111111));
    } else if(len == 3){
        unicode = ((unicode_t) ((uint8_t) s[0] & 0x0f)) << 12; // unicode = ((unicode_t) ((uint8_t) s[0] & 0b00001111)) << 12;
        unicode += ((unicode_t) ((uint8_t) s[1] & 0x3f)) << 6; // unicode += ((unicode_t) ((uint8_t) s[1] & 0b00111111)) << 6;
        unicode += ((unicode_t) ((uint8_t) s[2] & 0x3f)); // unicode += ((unicode_t) ((uint8_t) s[2] & 0b00111111));
    } else if(len == 4){
        unicode = ((unicode_t) ((uint8_t) s[0] & 0x07)) << 18; // unicode = ((unicode_t) ((uint8_t) s[0] & 0b00000111)) << 18;
        unicode += ((unicode_t) ((uint8_t) s[1] & 0x3f)) << 12; // unicode += ((unicode_t) ((uint8_t) s[1] & 0b00111111)) << 12;
        unicode += ((unicode_t) ((uint8_t) s[2] & 0x3f)) << 6; // unicode += ((unicode_t) ((uint8_t) s[2] & 0b00111111)) << 6;
        unicode += ((unicode_t) ((uint8_t) s[3] & 0x3f)); // unicode += ((unicode_t) ((uint8_t) s[3] & 0b00111111));
    }
    return unicode;
}

int32_t utf8_to_unicode_subset(const char * const s, int8_t * const len_p){
    unicode_t code = utf8_to_unicode(s, len_p);
    if(code == -1){return -1;}

   //  if(UNICODE_BASIC_LATIN_BEG <= code && code <= UNICODE_BASIC_LATIN_END){return UNICODE_BASIC_LATIN;}
    if(UNICODE_LATIN_UPPER_BEG <= code && code <= UNICODE_LATIN_UPPER_END){return UNICODE_LETTER;}
    if(UNICODE_LATIN_LOWER_BEG <= code && code <= UNICODE_LATIN_LOWER_END){return UNICODE_LETTER;}
    if(UNICODE_NUMBER_BEG <= code && code <= UNICODE_NUMBER_END){return UNICODE_NUMBER;}
    if(UNICODE_LATIN_PUNCTUATION1_BEG <= code && code <= UNICODE_LATIN_PUNCTUATION1_END){return UNICODE_PUNCTUATION;}
    if(UNICODE_LATIN_PUNCTUATION2_BEG <= code && code <= UNICODE_LATIN_PUNCTUATION2_END){return UNICODE_PUNCTUATION;}
    if(UNICODE_LATIN_PUNCTUATION3_BEG <= code && code <= UNICODE_LATIN_PUNCTUATION3_END){return UNICODE_PUNCTUATION;}
    if(UNICODE_LATIN_PUNCTUATION4_BEG <= code && code <= UNICODE_LATIN_PUNCTUATION4_END){return UNICODE_PUNCTUATION;}
    if(UNICODE_SUPPLEMENTAL_PUNCTUATION_BEG <= code && code <= UNICODE_SUPPLEMENTAL_PUNCTUATION_END){return UNICODE_PUNCTUATION;}
    if(UNICODE_CURRENCY_BEG <= code && code <= UNICODE_CURRENCY_END){return UNICODE_CURRENCY;}
    if(UNICODE_EMOTICON_BEG <= code && code <= UNICODE_EMOTICON_END){return UNICODE_EMOTICON;}
    if(UNICODE_PHONETIC_EXTENSION_BEG <= code && code <= UNICODE_PHONETIC_EXTENSION_END){return UNICODE_PHONETIC;}

    return UNICODE_OTHER;
}

char * utf8_normalise(char * const s){
    // int32_t count_letter = 0;
    int32_t count_number = 0;
    int32_t count_punctuation = 0;
    int32_t count_currency = 0;
    int32_t count_emoticon = 0;
    // int32_t count_phonetic = 0;
    // int32_t count_other = 0;

    int8_t found_hyphen = 0;
    int8_t found_non_hyphen = 0;

    int8_t len;

    int32_t index = 0;
    while(index < 64 && s[index] != '\0'){
        int32_t unicode_subset = utf8_to_unicode_subset(&s[index], &len);
        if(len == -1){index++; continue;}

        switch(unicode_subset){
//            case UNICODE_LETTER:
//                count_letter++;
//                break;
            case UNICODE_NUMBER:
                count_number++;
                break;
            case UNICODE_PUNCTUATION:
                count_punctuation++;
                if(len == 1){
                    if(s[index] == '-'){found_hyphen = 1;}
                    else {found_non_hyphen = 1;}
                }
                break;
            case UNICODE_CURRENCY:
                count_currency++;
                break;
            case UNICODE_EMOTICON:
                count_emoticon++;
                break;
//            case UNICODE_PHONETIC:
//                count_phonetic++;
//                break;
//            case UNICODE_OTHER:
//                count_other++;
//                break;
            default:
		break;
        }

        index += len;
    }

    if(found_hyphen && !found_non_hyphen){return s;}
    if(count_number >= 2){return "[NUMBER]";}
    if(count_emoticon >= 1){return "[EMOTICON]";}
    if(count_punctuation >= 2){return "[PUNCTUATION_GROUP]";}
    if(count_currency >= 1){return "[CURRENCY]";}

    return s;
}

