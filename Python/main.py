import zmq

import SerialHandler

context = zmq.Context()
zsock_t = context.socket(zmq.PULL)
zsock_t.bind("ipc:///tmp/tta")

sh = SerialHandler.SerialHandler('/dev/ttyACM0')
sh.startReadThread()

#sh.sendEncoderPublisher(False)
#sh.activateSensorPublisher(False)
sh.sendPidActivation(False)
sh.sendSafetyStopActivation(True)

print("Start listen")

global msg_code
global speed
global angle

msg_code = -1
speed = 0
angle = 0

def do_work():
    while True:
        global msg_code
        global speed
        global angle
        msg = zsock_t.recv()
        floats = msg.decode("utf-8").split(" ")
        msg_code = int(float(floats[0]))
        speed = float(floats[1])
        angle = float(floats[2])
        #print("Code: {}, Speed: {}, Angle: {}".format(msg_code,speed, angle))
        if msg_code == 0:
            sh.sendBrake(angle)
            return
        elif msg_code == 1:
            sh.sendMove(speed, angle)
        elif msg_code == 2:
            sh.sendBrake(angle)




do_work()
print('Shuting down')
sh.close()

