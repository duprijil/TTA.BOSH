#include <webcam.h>
#include <raylib.h>
#include <grab.h>
#include <unistd.h>

#include <zmqWrap.h>

#include <unistd.h>

#define MSG_EXIT 0
#define MSG_ACC 1
#define MSG_BRAKE 2

#define MAX_SPEED 26.0f // 28
#define MAX_ANGLE 24.0f // 26

#define SPEED 2.0f
#define ANGLE 2.0f

void send_info(zmqIPCPush* sock, float speed, float angle, char msg) {
    char buff[32];
    snprintf(buff, sizeof(buff), "%d %.2f %.2f",msg, speed, angle);
    sock->push_str(buff);
}

int main(void) {
	
    const RGBImage *frame_ptr = nullptr;
    Webcam webcam("/dev/video0", XRES, YRES);
    
    
    SetTraceLogLevel(LOG_NONE);
    InitWindow(XRES, YRES, "Camera monitor");
    SetTargetFPS(60);
    
    zmqIPCPush* zmqSocket = new zmqIPCPush("ipc:///tmp/tta");
    char c=0;
    float angle = 0.00f;
    float speed = 0.00f;
    int show_image = 1;
    int control = 1;
    while (!WindowShouldClose()) {
		if(IsKeyPressed(KEY_X)) {
			control *= -1;
		}
		if(control == 1)
		{
			if(IsKeyPressed(KEY_SPACE)) {
				speed = 0.00f;
				angle = 0.00f;
				c = MSG_BRAKE;
				send_info(zmqSocket, speed, angle, c);
			}
			else if(IsKeyDown(KEY_A)) {
				if(angle > -MAX_ANGLE) 
				{
					angle -= ANGLE;
					c = MSG_ACC;
					send_info(zmqSocket, speed, angle, c);	
				}
			}
			else if(IsKeyDown(KEY_D)) {
				if(angle < MAX_ANGLE)
				{
					angle += ANGLE;
					c = MSG_ACC;
					send_info(zmqSocket, speed, angle, c);
				}
			}        
			else if(IsKeyPressed(KEY_E)) {
				if(speed < MAX_SPEED)
				{
					speed += SPEED;
					c = MSG_ACC;
					send_info(zmqSocket, speed, angle, c);	
				}
			}
			else if(IsKeyPressed(KEY_Q)) {
				if(speed > -MAX_SPEED)
				{
					speed -= SPEED;
					c = MSG_ACC;
					send_info(zmqSocket, speed, angle,c);	
				}
			}
			else if(IsKeyPressed(KEY_W) || IsKeyPressed(KEY_S)) {
				angle = 0.00f;
				c = MSG_ACC;
				send_info(zmqSocket, speed, angle,c);
			}
			else if(IsKeyPressed(KEY_P)) {
				
				show_image *= -1;
			}
		}
		else
		{
			// Controlul automat
		}



		BeginDrawing();
		Texture2D texture;
		Image img;
		if(show_image == 1)
		{
			frame_ptr = &webcam.frame();
			Frame frame(*frame_ptr);
			img = LoadImageEx(frame.data,frame.w,frame.h);
			texture = LoadTextureFromImage(img);
			ClearBackground(BLACK);
			DrawTexture(texture,0,0,WHITE);
			DrawFPS(10,10);
		}
		EndDrawing();
		if(show_image == 1)
		{
			UnloadTexture(texture);
			UnloadImage(img);	
		}

    }

	speed = 0.00f;
	angle = 0.00f;
	c = MSG_EXIT;
	send_info(zmqSocket, speed, angle, c);

    CloseWindow();
    
    delete zmqSocket;
 
}
