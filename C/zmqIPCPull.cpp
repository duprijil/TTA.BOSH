#include <czmq.h>
#include "zmqWrap.h"

zmqIPCPull::zmqIPCPull(const char* path) {
	sock = zsock_new_pull(path);
}

char* zmqIPCPull::pull_str(void) {
	return zstr_recv(sock);
}

zmqIPCPull::~zmqIPCPull() {
	zsock_destroy(&sock);
}
