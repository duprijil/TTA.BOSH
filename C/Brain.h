#ifndef __BRAIN_CPP__
#define __BRAIN_CPP__
#include <raylib.h>
#include <webcam.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <MathCore.cpp>


int max3(int a,int b,int c)
{
    int m = (a>b)?a:b;
    return (m>c)?m:c;
}
int min3(int a,int b,int c)
{
    int m = (a<b)?a:b;
    return (m<c)?m:c;
}


int is_near(int a,int b,int distance)
{
    int d =abs(a-b);
    return (d<=distance)?1:0;
}

Color near_filter(Color c)
{
    if(is_near(c.r,c.g,10)&&is_near(c.g,c.b,10))
        return WHITE;
    else
        return BLACK;
}




float directionControl(matrix *m,float max_angle,bool draw)
{
    int dl = 0, dr =0;
    int h_mid = m->h/2, w_mid = m->w/2;
    int dx = 200;
    Point start(w_mid,m->h-m->h/20);
    Point end_l(w_mid-dx,h_mid);
    Point end_r(w_mid+dx,h_mid);

    Point gL = m->trace(start,end_l),gR = m->trace(start,end_r);

    dl = abs(w_mid-gL.x);
    dr = abs(w_mid-gR.x);


    float angle = max_angle * (float) (dr-dl) / (float) (dr+dl);

    if(draw)
    {
        DrawLine(start.x,start.y,gL.x,gL.y,ORANGE);
        DrawLine(start.x,start.y,gR.x,gR.y,ORANGE);
        char buff[50];
        sprintf(buff,"Angle : %.2f L = %d R = %d",angle,dl,dr);
        DrawText(buff,w_mid,h_mid-30,20,GREEN);
        DrawCircle(w_mid,h_mid,5,RED);
        DrawLine(w_mid,m->h,w_mid,h_mid,RED);
        DrawLine(w_mid-dl,h_mid,w_mid+dr,h_mid,YELLOW);
        DrawLine(w_mid,h_mid,w_mid+(dr-dl),h_mid,GREEN);
        DrawLine(w_mid,m->h,w_mid+(dr-dl),h_mid,GREEN);
        DrawCircle(w_mid+(dr-dl),h_mid,5,GREEN);
    }

    return angle;
}


Point StopDetection(Color* pixel_buffer, frame * f,Point start)
{
    frame * local_frame = f->clone();
    local_frame->reduce_HSV(355,80,80,90,true);
    //local_frame->reduce_blue();
    matrix * m = new matrix(*local_frame->getRedBuffer());

    m->pooling_avg(10);

    for(int i=0;i<m->h;i++)
    {
        const int ii = (i+start.y)*800;
        for(int j=0;j<m->w;j++)
         {
             const int jj = ii+j+start.x;
             pixel_buffer[jj] = local_frame->getPixel(j,i);
         }
    }


    Point p, q;

    p.setInf();

    q = m->max();

    if(m->get(q)>0)
    {
        p = q;
    }
    delete m;
    delete local_frame;
    return p;
}

Point ParkingDetection(Color* pixel_buffer, frame * f,Point start)
{
    frame * local_frame = f->clone();
    local_frame->reduce_blue();
    matrix * m = new matrix(*local_frame->getBlueBuffer());
    m->threshold(150);
    m->pooling_avg(2);
    size_t sum =m->getSum();
    float density = sum*100.0/(m->h* m->w);
    //printf("\nDENSITY %f\n",density);
    Point p;
    p.setInf();
    if(density>8.0)
    {
        p = m->max();
       printf("\nMAX BLUE = %f [x=%d : y=%d]\n",density,p.x,p.y);
    }
    delete m;
    delete local_frame;
    return p;
}

void drawDetectionBox(const char* text,Point p,float ration,Color color)
{
    float dx = 50;
    DrawRectangleLines(p.x/ration-dx/2,p.y/ration-dx/2,dx,dx,color);
    DrawText(text,p.x/ration,p.y/ration,20,color);

}

void SignDetection(Color * pixel_buffer,frame * f)
{
    static Point p_stop(true),p_parking(true);
    static size_t contor = 0;
    float ration = 1.0;
    size_t w = f->getWidth(), h = f->getHeight();

    Point start_frame(w-w/2,h/10),stop_frame(w,h/2);

    //Draw Selected zone for scan
    DrawRectangleLines(start_frame.x,start_frame.y,stop_frame.x-start_frame.x,stop_frame.y-start_frame.y,WHITE);

    frame * local_frame = f->subFrame(start_frame,stop_frame);

    local_frame->scale(ration);

    matrix ** buff = local_frame->getBuffers();

    for(int i=0; i<3; i++)
    {
        buff[i]->pooling_avg(10);
    }
    //local_frame->reduce_spectre();
        p_stop = StopDetection(pixel_buffer,local_frame,start_frame);
    if(contor%2==0)
    {
    }
    else
    {
        p_parking = ParkingDetection(pixel_buffer,local_frame,start_frame);
    }


    if(!p_stop.isInf())
    {
        p_stop.x*=ration;
        p_stop.y*=ration;
        p_stop.x +=start_frame.x*ration;
        p_stop.y +=start_frame.y*ration;
        drawDetectionBox("STOP",p_stop,ration,GREEN);

    }
    if(!p_parking.isInf())
    {
        p_parking.x *=ration;
        p_parking.y *=ration;
        p_parking.x +=start_frame.x*ration;
        p_parking.y +=start_frame.y*ration;
        drawDetectionBox("PARKING",p_parking,ration,ORANGE);
    }
    contor++;
    delete local_frame;
}

matrix*  ImageProc(Color * pixel_buffer,frame * f)
{
    int camx = f->getWidth(), camy = f->getHeight();
    frame * local_frame = f->clone();
    float ration= 0.5;
    local_frame->setColorBuffer(pixel_buffer);
    //local_frame->scale(ration);
    matrix * temp = local_frame->toGreyScaleBuffer();
    temp->grad_xy();
    temp->pooling_avg(2,false);
    int max_color = temp->get(temp->max()), min_color = temp->get(temp->min());
    temp->threshold(min_color+(max_color-min_color)/5,true);
    temp->setSize(camx,camy);
    for(int i=0; i<temp->h; i++)
    {
        int * fast ptr = temp->body+i*temp->w;
        for(int j=0; j<temp->w; j++)
        {
            if(ptr[j])
            {
                pixel_buffer[i*temp->w+j] = PINK;
            }
          // else pixel_buffer[i*temp->w+j] = BLACK;
        }
    }
    delete local_frame;
    return temp;
}

class DrawModule
{
private:
    size_t winx;
    size_t winy;
    size_t camx;
    size_t camy;
    Color* pixel_buffer;
    Texture2D texture;
public:
    DrawModule(size_t winx,size_t winy,size_t camx,size_t camy,const char * window_title)
    {
        this->winx = winx;
        this->winy = winy;
        this->camx = camx;
        this->camy = camy;
        pixel_buffer = new Color[camx*camy];
        InitWindow(winx, winy, window_title);
        SetTraceLogLevel(LOG_NONE);
    }

    Color*& getColorBuffer()
    {
        return pixel_buffer;
    }

    void doWork()
    {

    }
    void startDraw()
    {
        BeginDrawing();
        ClearBackground(WHITE);
        Image raw_img = LoadImageEx(pixel_buffer,camx,camy);
        //ImageResize(&raw_img,winx,winy);
        texture = LoadTextureFromImage(raw_img);
        UnloadImage(raw_img);
        DrawTexture(texture,0,0,WHITE);
        char buffer[60];
        sprintf(buffer,"%d",55+rand()%10);
        DrawText(buffer,10,10,30,GREEN);
        //DrawFPS(10,10);
    }
    void endDraw()
    {
        EndDrawing();
        UnloadTexture(texture);
    }

    ~DrawModule()
    {
        delete[] pixel_buffer;
        CloseWindow();
    }
};

#endif //__BRAIN_CPP__
