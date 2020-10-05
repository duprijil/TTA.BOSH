#include <czmq.h>
#include "zmqWrap.h"

zmqIPCPush::zmqIPCPush(const char* path) {
	sock = zsock_new_push(path);
}

void zmqIPCPush::push_str(const char* msg) {
	zstr_send(sock, msg);
}

zmqIPCPush::~zmqIPCPush() {
	zsock_destroy(&sock);
}
