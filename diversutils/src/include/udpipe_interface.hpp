#ifndef UDPIPE_INTERFACE_HPP
#define UDPIPE_INTERFACE_HPP

#include <iostream>
#include <string>
#include <sstream> // std::istringstream

#include <stdint.h>
#include <stdio.h>

#include "udpipe.h"
#include "udpipe_interface/size.h"

/*
std::string& ufal::udpipe::pipeline::get_input(){return input;}
std::string& ufal::udpipe::pipeline::get_tagger(){return tagger;}
std::string& ufal::udpipe::pipeline::get_parser(){return parser;}
std::string& ufal::udpipe::pipeline::get_output(){return output;}
*/

static ufal::udpipe::pipeline* global_pipeline;

extern "C" void ensure_proper_udpipe_pipeline_size(){
	if(UDPIPE_PIPELINE_SIZE != sizeof(ufal::udpipe::pipeline)){
		fprintf(stderr, "UDPIPE_PIPELINE_SIZE != sizeof(ufal::udpipe::pipeline); %i != %li\n", UDPIPE_PIPELINE_SIZE, sizeof(ufal::udpipe::pipeline));
		exit(1);
	}
}

/*
extern "C" void udpipe_pipeline_create(const char* model_name, const char* input, const char* tagger, const char* parser, const char* output, struct udpipe_pipeline* const pointer_pipeline){
	ufal::udpipe::pipeline local_pipeline(
		ufal::udpipe::model::load(model_name),
		std::string(input),
		std::string(tagger),
		std::string(parser),
		std::string(output)
	);

	std::cout << "input: " << std::string(input) << std::endl;
	std::cout << "tagger: " << std::string(tagger) << std::endl;
	std::cout << "parser: " << std::string(parser) << std::endl;
	std::cout << "output: " << std::string(output) << std::endl;

	std::cout << "get input:" << local_pipeline.get_input() << std::endl;

	memcpy(pointer_pipeline, std::addressof(local_pipeline), sizeof(ufal::udpipe::pipeline));
}
*/

extern "C" void udpipe_pipeline_create_global(const char* model_name, const char* input, const char* tagger, const char* parser, const char* output){
	std::cout << "Creating pipeline with '" << model_name << "'" << std::endl;
	static ufal::udpipe::pipeline local_pipeline = ufal::udpipe::pipeline(
		ufal::udpipe::model::load(model_name),
		std::string(input),
		std::string(tagger),
		std::string(parser),
		std::string(output)
	);
	global_pipeline = std::addressof(local_pipeline);
}

/*
extern "C" void udpipe_pipeline_print_global_info(){
	std::cout << "udpipe_pipeline_print_global_info:" << std::endl;
	std::cout << "get input:" << global_pipeline->get_input() << std::endl;
	std::cout << "get tagger:" << global_pipeline->get_tagger() << std::endl;
	std::cout << "get parser:" << global_pipeline->get_parser() << std::endl;
	std::cout << "get output:" << global_pipeline->get_output() << std::endl;
}
*/

extern "C" void udpipe_pipeline_process(const char* raw_txt, FILE** const pointer_file, void** const pointer_heap_char){
	std::string raw_txt_to_string(raw_txt);
	std::istringstream iss (raw_txt_to_string);
	// std::cout << "===================================================================\niss:\n" << iss.str() << std::endl;
	std::ostringstream oss("");
	std::string error_string;

	global_pipeline->process(iss, oss, error_string);

	if(error_string != ""){
		std::cerr << "error: " << error_string << std::endl;
		exit(1);
	}

	std::string str_from_oss = oss.str();
	size_t str_from_oss_length = str_from_oss.length();
	
	(*pointer_heap_char) = realloc((*pointer_heap_char), (str_from_oss_length + 1) * sizeof(char));
	if((*pointer_heap_char) == NULL){
		perror("realloc failed\n");
		exit(1);
	}
	memcpy((*pointer_heap_char), str_from_oss.c_str(), str_from_oss_length);
	((char*) (*pointer_heap_char))[str_from_oss_length] = '\0';
	(*pointer_file) = fmemopen((void*) (*pointer_heap_char), str_from_oss_length, "r");
}

#endif
