#ifndef __ZMQWRAP_H__
#define __ZMQWRAP_H__

#include <czmq.h>

class zmqIPCPush {
	private:
	zsock_t *sock;
	public:
	zmqIPCPush(const char* path);
	void push_str(const char* msg);
	~zmqIPCPush();
};

class zmqIPCPull {
	private:
	zsock_t *sock;
	public:
	zmqIPCPull(const char* path);
	char* pull_str(void);
	~zmqIPCPull();
};



#endif // __ZMQWRAP_H__
