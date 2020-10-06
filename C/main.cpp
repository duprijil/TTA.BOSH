#include <webcam.h>
#include <raylib.h>
#include <Brain.cpp>
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



class ImageProcessingModule
{
private:
    size_t w;
    size_t h;
    Webcam camera;
    frame * f;
    Color *data;
    zmqIPCPush* zmqSocket;
    DrawModule draw_m;
    bool source_select;
    char c=0;
    float angle = 0.00f;
    float speed = 0.00f;
    int show_image = 1;
    int control = 1;
public:
    ImageProcessingModule(size_t winx,size_t winy,size_t camx,size_t camy):camera("/dev/video0",camx,camy),
        draw_m(winx,winy,camx,camy,"Visual AI by Stratulat Stefan")
    {
        w = camx;
        h = camy;
        f = nullptr;
        source_select = true;
        camera.doPhotoAsync();
        zmqIPCPush* zmqSocket = new zmqIPCPush("ipc:///tmp/tta");
    }

    void send_info(zmqIPCPush* sock, float speed, float angle, char msg)
    {
        char buff[32];
        snprintf(buff, sizeof(buff), "%d %.2f %.2f",msg, speed, angle);
        sock->push_str(buff);
    }

    void KeyProcessing()
    {
        if(IsKeyPressed(KEY_X))
        {
            control *= -1;
        }
        if(control == 1)
        {
            if(IsKeyPressed(KEY_SPACE))
            {
                speed = 0.00f;
                angle = 0.00f;
                c = MSG_BRAKE;
                send_info(zmqSocket, speed, angle, c);
            }
            else if(IsKeyDown(KEY_A))
            {
                if(angle > -MAX_ANGLE)
                {
                    angle -= ANGLE;
                    c = MSG_ACC;
                    send_info(zmqSocket, speed, angle, c);
                }
            }
            else if(IsKeyDown(KEY_D))
            {
                if(angle < MAX_ANGLE)
                {
                    angle += ANGLE;
                    c = MSG_ACC;
                    send_info(zmqSocket, speed, angle, c);
                }
            }
            else if(IsKeyPressed(KEY_E))
            {
                if(speed < MAX_SPEED)
                {
                    speed += SPEED;
                    c = MSG_ACC;
                    send_info(zmqSocket, speed, angle, c);
                }
            }
            else if(IsKeyPressed(KEY_Q))
            {
                if(speed > -MAX_SPEED)
                {
                    speed -= SPEED;
                    c = MSG_ACC;
                    send_info(zmqSocket, speed, angle,c);
                }
            }
            else if(IsKeyPressed(KEY_W) || IsKeyPressed(KEY_S))
            {
                angle = 0.00f;
                c = MSG_ACC;
                send_info(zmqSocket, speed, angle,c);
            }
            else if(IsKeyPressed(KEY_P))
            {

                show_image *= -1;
            }
        }
        else
        {
            // Controlul automat
        }
    }
    void updateFrame()
    {
        if(f!=nullptr)
            delete f;
        if(source_select)
        {


            //while(!camera.isDone());
            f = new frame(camera.frame());
            //camera.doPhotoAsync();
        }
        else
        {
            f = new frame(data,w,h);
        }
    }
    void doWork()
    {
        SetTraceLogLevel(LOG_NONE);
        SetTargetFPS(60);
        size_t contor = 0;
        char buff[256];
        float img_proc_time =0,sign_proc_time=0,dir_proc_time=0;
        while (!WindowShouldClose())
        {
            KeyProcessing();
            updateFrame();

            f->normalize();

            Clock bench;

            matrix * temp = ImageProc(draw_m.getColorBuffer(),f);

            img_proc_time = bench.getElapsed();

            bench.update();
            draw_m.startDraw();
            //SignDetection(draw_m.getColorBuffer(),f);


            sign_proc_time = bench.getElapsed();
            bench.update();

            directionControl(temp,30,true);

            dir_proc_time = bench.getElapsed();

            sprintf(buff,"IMG_PROC = %.3fs\nSIGN_PROC = %.3fs\nDIR_PROC = %.3fs\n",img_proc_time,sign_proc_time,dir_proc_time);

            DrawText(buff,10,50,30,GREEN);

            draw_m.endDraw();
            delete temp;
            contor++;
        }
    }

    ~ImageProcessingModule()
    {
        speed = 0.00f;
        angle = 0.00f;
        c = MSG_EXIT;
        send_info(zmqSocket, speed, angle, c);
        delete zmqSocket;
        delete[] data;
        if(f!=nullptr)
            delete f;
    }

};




int main(void)
{
    int winx = 800, winy = 600, camx= 800, camy = 600;

    ImageProcessingModule img_m(winx,winy,camx,camy);
    img_m.doWork();
}
