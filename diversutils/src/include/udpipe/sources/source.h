#ifndef UDPIPE_SOURCE_H
#define UDPIPE_SOURCE_H

struct udpipe_source {
    const char * const url_repository;
    const char * const url_download;
    const char isolang[4];
    const char md5hex[33];
    const char * const filename;
};

#endif
