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

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef PCRE2_CODE_UNIT_WIDTH
#define PCRE2_CODE_UNIT_WIDTH 0
#include <pcre2.h>
#endif

#include "filter.h"
#include "logging.h"

int32_t filter_compile(struct filter * const f, const uint8_t quiet){
    int errorcode;
    PCRE2_SIZE erroroffset;
    const size_t log_bfr_size = 1024;
    char log_bfr[log_bfr_size];
    
    if(!quiet){
        memset(log_bfr, '\0', log_bfr_size);
        snprintf(log_bfr, log_bfr_size, "Attempting to compile: %s <-> %s", f->replacement, f->pattern_match);
        info_format(__FILE__, __func__, __LINE__, log_bfr);
    }

    f->length = strlen((char*) f->pattern_match);
    f->rlength = strlen((char*) f->replacement);
    f->regex = pcre2_compile(
        f->pattern_match,
        PCRE2_ZERO_TERMINATED,
        f->options_compile,
        &errorcode,
        &erroroffset,
        NULL
    );
    if(f->regex == NULL){
        const int32_t bfr_size = 256;
        unsigned char bfr[bfr_size];
        memset(bfr, '\0', bfr_size);
        pcre2_get_error_message(errorcode, bfr, bfr_size);
        fprintf(stderr, "error: %s\n", bfr);
        return 1;
    }
    return 0;
}

void filter_free(struct filter * const f){
    pcre2_code_free(f->regex);
}

struct filter filters[] = {
    #if ENABLE_FILTER_XML == 1
    { .regex = NULL, .pattern_match = (unsigned char*) "</?[a-z0-9]{1,32}+(\\s+[a-z0-9]{1,32}+=\"[^\"]{0,128}+\"){0,16}+\\s?/?>", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[XML]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    #endif
    #if ENABLE_FILTER_PATH == 1
    { .regex = NULL, .pattern_match = (unsigned char*) "(?<!/)(/[a-z0-9\\-\\._]{1,64}+){1,32}+/?", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[PATH]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    #endif
    #if ENABLE_FILTER_EMAIL == 1
    { .regex = NULL, .pattern_match = (unsigned char*) "[a-z0-9\\-_\\.]{1,32}+@([a-z0-9\\-_]{1,32}+\\.){1,4}+[a-z0-9]{1,6}+", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[EMAIL]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    #endif
    #if ENABLE_FILTER_URL == 1
    { .regex = NULL, .pattern_match = (unsigned char*) "\n\
        ([a-z0-9]{1,6}://)?\n\
        (\n\
            \\[EMAIL\\]\n\
            |\n\
            ([a-z0-9\\-_\\.]{1,32}@)?\n\
            ([a-z0-9\\-_]{1,16}\\.){1,4}[a-z0-9]{1,6}\n\
        )\n\
        (:\\d{,5})?\n\
        (\\[PATH\\]|(/[^/\\&;\\s]{1,64}+){1,16}+)?\n\
        (\\?([^=]{1,32}+=[^\\&;]{1,32}+[\\&;]?){,32}+)?\n\
        (\\#\\S{,64}+)?\n\
        ", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS | PCRE2_EXTENDED, .replacement = (unsigned char*) "[URL]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    /*
    { .regex = NULL, .pattern_match = (unsigned char*) "\n\
    # scheme (see https://www.iana.org/assignments/uri-schemes/uri-schemes.xhtml) \n\
    ((aaa|aaas|about|acap|acct|acd|acr|adiumxtra|adt|afp|afs|aim|amss|android|appdata|apt|ar|ark|at|attachment|aw|barion|bb|beshare|bitcoin|bitcoincash|blob|bolo|brid|browserext|cabal|calculator|callto|cap|cast|casts|chrome|chrome-extension|cid|coap|coap\\+tcp|coap\\+ws|coaps|coaps\\+tcp|coaps\\+ws|com-eventbrite-attendee|content|content-type|crid|cstr|cvs|dab|dat|data|dav|dhttp|diaspora|dict|did|dis|dlna-playcontainer|dlna-playsingle|dns|dntp|doi|dpp|drm|drop|dtmi|dtn|dvb|dvx|dweb|ed2k|eid|elsi|embedded|ens|ethereum|example|facetime|fax|feed|feedready|fido|file|filesystem|finger|first-run-pen-experience|fish|fm|ftp|fuchsia-pkg|geo|gg|git|gitoid|gizmoproject|go|gopher|graph|grd|gtalk|h323|ham|hcap|hcp|hs20|http|https|hxxp|hxxps|hydrazone|hyper|iax|icap|icon|im|imap|info|iotdisco|ipfs|ipn|ipns|ipp|ipps|irc|irc6|ircs|iris|iris\\.beep|iris\\.lwz|iris\\.xpc|iris\\.xpcs|isostore|itms|jabber|jar|jms|keyparc|lastfm|lbry|ldap|ldaps|leaptofrogans|lid|lorawan|lpa|lvlt|machineProvisioningProgressReporter|magnet|mailserver|mailto|maps|market|matrix|message|microsoft\\.windows\\.camera|microsoft\\.windows\\.camera\\.multipicker|microsoft\\.windows\\.camera\\.picker|mid|mms|modem|mongodb|moz|ms-access|ms-appinstaller|ms-browser-extension|ms-calculator|ms-drive-to|ms-enrollment|ms-excel|ms-eyecontrolspeech|ms-gamebarservices|ms-gamingoverlay|ms-getoffice|ms-help|ms-infopath|ms-inputapp|ms-launchremotedesktop|ms-lockscreencomponent-config|ms-media-stream-id|ms-meetnow|ms-mixedrealitycapture|ms-mobileplans|ms-newsandinterests|ms-officeapp|ms-people|ms-powerpoint|ms-project|ms-publisher|ms-recall|ms-remotedesktop|ms-remotedesktop-launch|ms-restoretabcompanion|ms-screenclip|ms-screensketch|ms-search|ms-search-repair|ms-secondary-screen-controller|ms-secondary-screen-setup|ms-settings|ms-settings-airplanemode|ms-settings-bluetooth|ms-settings-camera|ms-settings-cellular|ms-settings-cloudstorage|ms-settings-connectabledevices|ms-settings-displays-topology|ms-settings-emailandaccounts|ms-settings-language|ms-settings-location|ms-settings-lock|ms-settings-nfctransactions|ms-settings-notifications|ms-settings-power|ms-settings-privacy|ms-settings-proximity|ms-settings-screenrotation|ms-settings-wifi|ms-settings-workplace|ms-spd|ms-stickers|ms-sttoverlay|ms-transit-to|ms-useractivityset|ms-virtualtouchpad|ms-visio|ms-walk-to|ms-whiteboard|ms-whiteboard-cmd|ms-word|msnim|msrp|msrps|mss|mt|mtqp|mumble|mupdate|mvn|mvrp|mvrps|news|nfs|ni|nih|nntp|notes|num|ocf|oid|onenote|onenote-cmd|opaquelocktoken|openid|openpgp4fpr|otpauth|p1|pack|palm|paparazzi|payment|payto|pkcs11|platform|pop|pres|prospero|proxy|psyc|pttp|pwid|qb|query|quic-transport|redis|rediss|reload|res|resource|rmi|rsync|rtmfp|rtmp|rtsp|rtsps|rtspu|sarif|secondlife|secret-token|service|session|sftp|sgn|shc|shttp\\ \\(OBSOLETE\\)|sieve|simpleledger|simplex|sip|sips|skype|smb|smp|sms|smtp|snews|snmp|soap\\.beep|soap\\.beeps|soldat|spiffe|spotify|ssb|ssh|starknet|steam|stun|stuns|submit|svn|swh|swid|swidpath|tag|taler|teamspeak|tel|teliaeid|telnet|tftp|things|thismessage|thzp|tip|tn3270|tool|turn|turns|tv|udp|unreal|upt|urn|ut2004|uuid-in-package|v-event|vemmi|ventrilo|ves|videotex|view-source|vnc|vscode|vscode-insiders|vsls|w3|wais|wcr|web\\+ap|web3|webcal|wifi|wpid|ws|wss|wtai|wyciwyg|xcon|xcon-userid|xfire|xftp|xmlrpc\\.beep|xmlrpc\\.beeps|xmpp|xrcp|xri|ymsgr|z39\\.50|z39\\.50r|z39\\.50s):)?\n\
    # authority \n\
    (//\n\
        # userinfo \n\
        #(([a-zA-Z0-9]([a-zA-Z0-9\\-_\\.]*[a-zA-Z0-9])?)@)?\n\
        ([^@\\s]+@)?\n\
        # host \n\
        (\n\
            # ipv4 \n\
            (\\d{1,3}(\\.\\d{1,3}){3}) |\n\
            # ipv6 \n\
            (([0-9a-fA-F]{4}:){1,7}([0-9a-fA-F]{4}|:)) |\n\
            #([0-9a-fA-F]{,4}:){2,8} |\n\
            # domain \n\
                # subdomains and main part \n\
                #([a-zA-Z0-9\\-_]+\\.)+\n\
                ([^\\.]+\\.)+\n\
                # domain type \n\
                [a-zA-Z0-9]+\n\
        )\n\
        # port \n\
        (:\\d+)?\n\
    )?\n\
    # path \n\
    #(/[a-zA-Z0-9\\-_\\.]+)*\n\
    (/[^/\\#\\?]+)*\n\
    # query \n\
    (\\?\n\
        (\n\
            # key \n\
            #([a-zA-Z%0-9]+)\n\
            [^=]+\n\
            =\n\
            # value \n\
            #([a-zA-Z%0-9]+)\n\
            [^\\&;]+\n\
            # and \n\
            [\\&;]?\n\
        )*\n\
    )?\n\
    # fragment \n\
    (\\#\n\
        [a-zA-Z%0-9]*\n\
    )?\n\
    ", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_EXTENDED, .replacement = (unsigned char*) "[URL]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    */
    #endif
    #if ENABLE_FILTER_ALPHANUM == 1
    // { .regex = NULL, .pattern_match = (unsigned char*) "[a-z0-9\\+\\-\\.]{,32}[0-9][a-z0-9\\+\\-\\.]{,32}", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[ALPHANUM]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    // { .regex = NULL, .pattern_match = (unsigned char*) "\\w{,32}[0-9]\\w{,32}", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[ALPHANUM]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    { .regex = NULL, .pattern_match = (unsigned char*) "\\w*?[0-9]\\w*+", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[ALPHANUM]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    #endif
    #if ENABLE_FILTER_LONG == 1
    { .regex = NULL, .pattern_match = (unsigned char*) "\\S{30,}+", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL, .replacement = (unsigned char*) "[LONG]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    #endif
    #if ENABLE_FILTER_NON_FRENCH == 1
    { .regex = NULL, .pattern_match = (unsigned char*) "\\w{0,8}[^0-9a-zàâäéèêëìîïòôöùûüç\\s&~\"#'\\{\\(\\[\\|\\|\\`_\\\\\\^@\\)\\]°\\+=\\}\\$£%µ\\*§!:/;\\.,\\?<>]\\w++", .options_compile = PCRE2_UCP | PCRE2_UTF | PCRE2_DOTALL | PCRE2_CASELESS, .replacement = (unsigned char*) "[NON_FRENCH]", .options_replace = PCRE2_SUBSTITUTE_LITERAL | PCRE2_SUBSTITUTE_GLOBAL | PCRE2_SUBSTITUTE_OVERFLOW_LENGTH, },
    #endif
};

const size_t filter_cardinality = sizeof(filters) / sizeof(struct filter);

int32_t filter_ready(const uint8_t quiet){
    size_t i;
    for(i = 0 ; i < filter_cardinality ; i++){
        if(filter_compile(&(filters[i]), quiet) != 0){
            fprintf(stderr, "Failed to call filter_compile for %s\n", filters[i].pattern_match);
            return 1;
        }
    }
    return 0;
}

void filter_release(void){
    size_t i;
    for(i = 0 ; i < filter_cardinality ; i++){
        filter_free(&(filters[i]));
    }
}

int32_t filter_substitute_all(char ** const s, size_t * const s_len){
    size_t i;
    PCRE2_SPTR subject = (unsigned char*) (*s);
    // PCRE2_SIZE length = strlen(*s);
    pcre2_match_data * match_data = NULL;
    pcre2_match_context * mcontext = NULL;
    // PCRE2_SIZE outputlength = length + 1;
    PCRE2_SIZE outputlength = (*s_len) + 1;
    PCRE2_UCHAR * outputbuffer = malloc(outputlength);
    if(outputbuffer == NULL){
        perror("alloc failed\n");
        return 1;
    }
    // memset(outputbuffer, '\0', outputlength); // ?

    int num_sub = 0;

    for(i = 0 ; i < filter_cardinality ; i++){
        assert(outputbuffer != NULL);

        PCRE2_SIZE startoffset = 0;

        substitution:

        num_sub = pcre2_substitute(
            filters[i].regex,
            subject,
            PCRE2_ZERO_TERMINATED,
            startoffset,
            filters[i].options_replace,
            match_data,
            mcontext,
            filters[i].replacement,
            PCRE2_ZERO_TERMINATED,
            outputbuffer,
            &outputlength
        );
        if(num_sub == PCRE2_ERROR_NOMEMORY){
            outputlength++;
            void * realloc_ptr = realloc(outputbuffer, outputlength);
            if(realloc_ptr == NULL){
                perror("alloc failed\n");
                return 1;
            }
            outputbuffer = realloc_ptr;
            goto substitution;
        }

        if(num_sub < 0){
            fprintf(stderr, "failed to call pcre2_substitute for %s; return: %i\n", filters[i].pattern_match, num_sub);
            /*
            memset(outputbuffer, '\0', outputlength);
            pcre2_get_error_message(num_sub, outputbuffer, outputlength);
            fprintf(stderr, "error message: %s\n", outputbuffer);
            */
            #ifndef NDEBUG
            const size_t n = 1024;
            PCRE2_UCHAR z[n];
            memset(z, '\0', n * sizeof(PCRE2_UCHAR));
            pcre2_get_error_message(num_sub, z, n);
            fprintf(stderr, "error message: %s\n", z);
            #endif
            return 1;
        }
        if(num_sub != 0){
            // DO NOT REMOVE - frees *s, moves outputbuffer into *s, reallocates outputbuffer
            free(*s);
            (*s) = (char*) outputbuffer;
            #if PCRE2_CODE_UNIT_WIDTH == 0
            (*s_len) = strlen(outputbuffer);
            #else
            (*s_len) = outputlength * (PCRE2_CODE_UNIT_WIDTH >> 3);
            #endif
            outputbuffer = NULL;
            if(i < filter_cardinality - 1){
                subject = (PCRE2_SPTR) (*s);
                
                // length = outputlength;
                if(outputlength == 0){break;}
                // outputlength = length + 1;
                // outputbuffer = malloc(outputlength);
                outputbuffer = malloc(outputlength + 1);
                if(outputbuffer == NULL){
                    perror("alloc failed\n");
                    return 1;
                }
                // memset(outputbuffer, '\0', outputlength); // ?
            }
            

            /* // DO NOT REMOVE - reallocs *s and keeps outputbuffer
            const size_t realloc_size = (outputlength * (PCRE2_CODE_UNIT_WIDTH >> 3)) + 1;
            (*s) = realloc(*s, realloc_size);
            if(*s == NULL){
                perror("realloc failed\n");
                return 1;
            }
            memcpy(*s, outputbuffer, realloc_size - 1);
            (*s)[realloc_size - 1] = '\0';

            subject = (unsigned char*) *s;
            */
        } /* else {
            memset(outputbuffer, '\0', outputlength);
        }*/
    }

    if(outputbuffer != NULL){
        free(outputbuffer);
    }

    return 0;
}
