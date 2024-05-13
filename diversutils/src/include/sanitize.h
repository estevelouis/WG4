#ifndef SANITIZE_H
#define SANITIZE_H

int32_t sanitize_for_shell(const char* const s, const size_t n, char** const z){
	size_t alloc_size = n + 1024;
	int32_t i, j;

	(*z) = malloc(alloc_size * sizeof(char));
	if((*z) == NULL){
		perror("malloc failed\n");
		return 1;
	}
	memset((*z), '\0', alloc_size * sizeof(char));

	j = 0;
	for(i = 0 ; i < n ; i++){
		if(j + 2 >= alloc_size){alloc_size += 1024; (*z) = realloc((*z), alloc_size * sizeof(char)); if((*z) == NULL){perror("failed to realloc\n"); return 1;}}
		switch(s[i]){
			case '\\':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = '\\';
				j++;
				break;
			case '\n':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = 'n';
				j++;
				break;
			case '\r':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = 'r';
				j++;
				break;
			case '"':
				(*z)[j] = '\\';
				j++;
				(*z)[j] = '"';
				j++;
				break;
			default:
				(*z)[j] = s[i];
				j++;
				break;
		}
	}
	return 0;
}

#endif
