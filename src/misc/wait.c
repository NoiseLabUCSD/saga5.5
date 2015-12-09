#include <sys/types.h>
#include <sys/wait.h>
int wait_(status) {
	return wait(&status);
}

