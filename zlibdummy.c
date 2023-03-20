#include <stdio.h>

FILE * gzopen(const char *path, const char *mode) {
  char mode1[2];
  // the T flag means do not compress; it should be ignored
  // by fopen, but we never know ...
  mode1[0] = mode[0];
  if( mode[1] == 'T') mode1[1]='b';
  else mode1[1]=mode[1];
  
  return fopen(path,mode1);
}

int gzclose(FILE *fp) {
  return fclose(fp);
}

char * gzgets(FILE *stream, char *s, int size) {
  return fgets(s,size,stream);
}

size_t gzread(FILE *stream, void *ptr, size_t size) {
  return fread(ptr, (size_t) 1, size, stream);
}

int gzseek(FILE *stream, long offset, int whence) {
  return fseek(stream, offset, whence);
}

size_t gzwrite(FILE *stream,const void *ptr,size_t size) {
  return fwrite(ptr,(size_t) 1,size,stream);
}

int gzputc(FILE *stream, int c) {
  return fputc(c, stream);
}


int gzrewind(FILE *stream) {
  rewind(stream);
  return 0;
}
