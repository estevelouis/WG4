#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "jsonl/parser.h"

#ifndef JSONL_CONTENT_KEY
#define JSONL_CONTENT_KEY "text"
#endif

#ifndef TARGET_NUMBER_DOCUMENTS
#define TARGET_NUMBER_DOCUMENTS 1000000
#endif

int32_t main(int32_t argc, char** argv){
	int32_t argi;
	struct jsonl_document_iterator jdi;
	uint64_t num_documents;
	FILE* f_out;

	if(argc < 3){
		perror("error: argc < 3\n");
		return 1;
	}

	f_out = fopen(argv[1], "w");
	if(f_out == NULL){
		fprintf(stderr, "failed to open %s for write\n", argv[1]);
		return 1;
	}

	num_documents = 0;

	for(argi = 2 ; argi < argc ; argi++){
		printf("%i: %s\n", argi, argv[argi]);
		jdi = (struct jsonl_document_iterator) {};
		if(create_jsonl_document_iterator(&jdi, argv[argi], JSONL_CONTENT_KEY) != 0){
			perror("failed to call create_jsonl_document_iterator\n");
			return 1;
		}

		if(iterate_jsonl_document_iterator(&jdi) != 0){
			perror("failed to call iterate_jsonl_document_iterator\n");
			free_jsonl_document_iterator(&jdi);
			return 1;
		}

		while(!(jdi.file_is_done)){
			fprintf(f_out, "%s\n", jdi.current_document.text);

			num_documents++;
	
			if(num_documents < UINT64_MAX){
				if(num_documents >= TARGET_NUMBER_DOCUMENTS){
					printf("reached num_documents >= TARGET_NUMBER_DOCUMENTS\n");
					break;
				}
			} else {
				printf("reached num_documents == UINT64_MAX\n");
				break;
			}

			if(iterate_jsonl_document_iterator(&jdi) != 0){
				perror("failed to call iterate_jsonl_document_iterator\n");
				free_jsonl_document_iterator(&jdi);
				return 1;
			}
		}

		free_jsonl_document_iterator(&jdi);

		if(num_documents >= TARGET_NUMBER_DOCUMENTS){break;}
	}

	fclose(f_out);

	printf("number of effectively processed documents: %i\n", num_documents);

	return 0;
}
