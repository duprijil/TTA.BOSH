

#include <cmath>
#include <raylib.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <utility>
#include <string.h>
#include <RGBImage.h>
#define RAD(A)  (M_PI*((double)(A))/180.0)
typedef unsigned char uint8_t;

class Clock
{
private:
    float current_time;

    float CPS()
    {
        return (float)clock()/(float)CLOCKS_PER_SEC;
    }

public:
    Clock()
    {
        current_time = CPS();
    }
    bool isElapsed(float seconds)
    {
        float dif_time = CPS();
        dif_time -=current_time;
        if(dif_time>seconds)
            return true;
        else
            return false;
    }

    bool onElapsed(float seconds)
    {
        if(isElapsed(seconds))
        {
            update();
            return true;
        }
        else
            return false;
    }

    float getElapsed()
    {
        return CPS()-current_time;
    }

    void update()
    {
        current_time = CPS();
    }
};


struct Point
{
    int x;
    int y;
    Point(int x = 0,int y = 0)
    {
        this->x = x;
        this->y = y;
    }
    Point(bool INF)
    {
        if(INF)
        {
            x = INT_MAX;
            y = INT_MAX;
        }
    }

    void replaceOnX(Point& p)
    {
        if(x>p.x)
        {
            int temp = x;
            x = p.x ;
            p.x = temp;
            temp= y;
            y=p.y;
            p.y=temp;
        }
    }

    void replaceOnY(Point& p)
    {
        if(y>p.y)
        {
            int temp_y = y;
            y = p.y ;
            p.y = temp_y;
        }
    }

    void setInf()
    {
        x = INT_MAX;
        y = INT_MAX;
    }

    bool isInf()
    {
        if((x == INT_MAX) && (y==INT_MAX))
            return true;
        else
            return false;
    }
    double distance(Point &p)
    {
        double dx = p.x-x;
        double dy = p.y-y;
        return sqrt(dx*dx+dy*dy);
    }
};

struct Line
{
    Point start;
    Point stop;
    Line()
    {
        start.x = 0;
        start.y = 0;
        stop.x = 0;
        stop.y = 0;
    }
    Line(Point &start,Point &stop)
    {
        this->start = start;
        this->stop = stop;
    }

    float getLength()
    {
        return start.distance(stop);
    }

};

#define fast const __restrict

struct matrix
{
    int * __restrict body;
    int w;
    int h;
    matrix()
    {
        body = nullptr;
        w =0;
        h =0;
    }
    matrix(int w,int h)
    {
        this->w = w;
        this->h = h;
        body = new int[w*h];
    }

    matrix(const matrix &m):matrix(m.w,m.h)
    {
        memcpy(body,m.body,w*h*sizeof(int));
    }

    matrix(const uint8_t * pixel_buffer,int chanel,int w,int h):matrix(w,h)
    {
        int k=0;
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
            {
                ptr[j] = pixel_buffer[k+chanel];
                k+=3;
            }
        }
    }

    ~matrix()
    {
        delete[] body;
    }

    void set(int value)
    {
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
            {
                ptr[j] = value;
            }
        }
    }

    void setPixelBuffer(uint8_t *pixel_buff)
    {
        int k = 0 ;
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
            {
                register const int value = ptr[j];
                pixel_buff[k] = value;
                pixel_buff[k+1] = value;
                pixel_buff[k+2] = value;
                k+=3;
            }
        }
    }

    Point trace(Point start,Point end)
    {
        Point local(end.x,end.y);
        double panta = ((double)(start.y-end.y)/(double)(start.x-end.x));
        if(start.x<end.x)
        {
            for(int x = start.x; x<end.x; x++)
            {
                int y = panta*(x-start.x) + start.y;

                int * fast ptr = body+y*w;

                if((x>=0)&&(y>=0)&&(x<w)&&(y<h))
                    if(ptr[x]!=0)
                    {
                        local.x=x;
                        local.y=y;
                        return local;
                    }
            }
        }
        else
        {
            for(int x = start.x; x>end.x; x--)
            {
                int y = panta*(x-start.x) + start.y;

                int * fast ptr = body+y*w;

                if((x>=0)&&(y>=0)&&(x<w)&&(y<h))
                    if(ptr[x]!=0)
                    {
                        local.x=x;
                        local.y=y;
                        DrawLine(start.x,start.y,local.x,local.y,RED);
                        return local;
                    }
            }
        }
        return local;
    }

    void setSize(int w,int h)
    {
        if((w>0)&&(h>0)&&((w!=this->w)||(h!=this->h)))
        {
            int * const new_mat = new int[w*h];

            const int old_w = this->w, old_h = this->h;

            register const float step_x=(float)old_w/(float)w;
            register const float step_y=(float)old_h/(float)h;
            float x=0.0f,y=0.0f;

            for(int i=0; i<h; i++)
            {
                x=0.0f;
                int * fast dest = new_mat+i*w;
                int * fast src = body+((int)trunc(y))*w;
                for(int j=0; j<w; j++)
                {
                    dest[j]=src[(int)trunc(x)];
                    x+=step_x;
                }
                y+=step_y;
            }
            delete[] body;
            body = new_mat;
            this->w = w;
            this->h = h;
        }
    }

    void scale(float ration)
    {
        setSize(w*ration,h*ration);
    }

    matrix * clone()
    {
        matrix * m = new matrix(*this);
        return m;
    }
    matrix * clone(matrix &m)
    {
        matrix *local = new matrix(m);
        return local;
    }

    void grad_x()
    {
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<(w-1); j++)
            {
                ptr[j]=abs(ptr[j]-ptr[j+1]);
            }
        }
        for(int i=0; i<h; i++)
            body[i*w+w-1] = 0;

    }
    void grad_y()
    {
        for(int i=0; i<(h-1); i++)
        {
            int * fast ptr = body+i*w;
            int * fast ptr_down = body+(i+1)*w;
            for(int j=0; j<w; j++)
            {
                ptr[j]=abs(ptr[j]-ptr_down[j]);
            }
        }
        const int index = (h-1)*w;
        for(int i=0; i<w; i++)
            body[index+i] = 0;

    }

    void grad_xy()
    {
        for(int i=0; i<(h-1); i++)
        {
            int * fast ptr = body+i*w;
            int * fast ptr_down = body+(i+1)*w;
            for(int j=0; j<(w-1); j++)
            {
                ptr[j]=abs(ptr[j]-ptr[j+1])+abs(ptr[j]-ptr_down[j]);
            }
        }
        for(int i=0; i<h; i++)
            body[i*w+w-1] = 0;

        const int index = (h-1)*w;
        for(int i=0; i<w; i++)
            body[index+i] = 0;

    }

    void setColorBuffer(Color * pixel_buffer)
    {
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
            {
                const int index = i*w+j, value = ptr[j];
                pixel_buffer[index].r = value;
                pixel_buffer[index].g = value;
                pixel_buffer[index].b = value;
                pixel_buffer[index].a = 0xFF;
            }
        }
    }

    void pooling_avg(int level,bool upScale=true)
    {
        if((level<w)&&(level<h)&&(level>1))
        {
            int old_w = w, old_h = h;
            setSize(w-w%level,h-h%level);
            int sum,i,j;
            for(i=0; i<(h-h%level); i+=level)
            {
                for(j=0; j<(w-w%level); j+=level)
                {
                    sum=0;
                    for(int k=i; k<(i+level); k++)
                    {
                        int * fast ptr = body+k*w;
                        for(int m=j; m<(j+level); m++)
                            sum+=ptr[m];
                    }
                    sum=sum/(level*level);
                    for(int k=i; k<(i+level); k++)
                    {
                        int * fast ptr = body+k*w;
                        for(int m=j; m<(j+level); m++)
                            ptr[m] = sum;
                    }
                }
            }
            if(upScale)
                setSize(old_w,old_h);
        }
    }
    void pooling_max(int level,bool upScale=true)
    {
        if((level<w)&&(level<h)&&(level>1))
        {
            int old_w = w, old_h = h;
            setSize(w-w%level,h-h%level);
            int maximum=0,i,j;
            for(i=0; i<(h-h%level); i+=level)
            {
                for(j=0; j<(w-w%level); j+=level)
                {
                    maximum =0;
                    for(int k=i; k<(i+level); k++)
                    {
                        int * fast ptr = body+k*w;
                        for(int m=j; m<(j+level); m++)
                        {
                            if(maximum<ptr[m])
                                maximum=ptr[m];
                        }
                    }
                    for(int k=i; k<(i+level); k++)
                    {
                        int * fast ptr = body+k*w;
                        for(int m=j; m<(j+level); m++)
                        {
                            ptr[m]= maximum;
                        }
                    }
                }
            }
            if(upScale)
                setSize(old_w,old_h);
        }
    }


    unsigned int getSum()
    {
        unsigned int sum = 0;
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
                sum+=ptr[j];
        }
        return sum;
    }

    Point min()
    {
        int min_val = body[0];

        Point p;

        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
                if(min_val >ptr[j])
                {
                    min_val = ptr[j];
                    p.x=j;
                    p.y=i;
                }
        }

        return p;
    }

    int get(Point p)
    {
        return body[p.y*w+p.x];
    }

    int get(size_t i,size_t j)
    {
        return body[i*w+j];
    }
    Point max()
    {
        Point p;
        int max_val = body[0];

        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
                if(max_val < ptr[j])
                {
                    max_val = ptr[j];
                    p.x = j;
                    p.y = i;
                }
        }
        return p;
    }

    void set(matrix &m)
    {
        h=m.h;
        w=m.w;
        delete[] body;
        body = new int[w*h];
        memcpy(body,m.body,w*h*sizeof(int));
    }

    matrix* conv(matrix *kernel,bool as_float_kernel, bool padding = true)
    {
        if(kernel->h!=kernel->w)
        {
            return nullptr;
        }
        if((kernel->h%2==0))
        {
            return nullptr;
        }
        int dx = (padding)?kernel->h/2:0;
        int dy = (padding)?kernel->h/2:0;
        int new_w = (padding)?w:w-kernel->h;
        int new_h = (padding)?h:h-kernel->h;
        int div = kernel->getSum();
        if(div==0)
            div=1;

        matrix * local = new matrix(new_w,new_h);

        for(int i=0; i<new_h; i++)
        {
            for(int j=0; j<new_w; j++)
            {
                int c_var = 0;
                int ii = i - dy;

                int * fast ptr = body+ii*w;

                if(as_float_kernel)
                {
                    float acum = 0.0f;
                    for(int n=0; n<kernel->h; n++)
                    {
                        const int * fast kernel_body_n = kernel->body+n*kernel->w;
                        int jj = j - dx;

                        for(int m=0; m<kernel->h; m++)
                        {
                            const int kernel_body_nm = kernel_body_n[m];
                            if(kernel_body_nm)
                            {
                                if((ii>=0)&&(ii<h)&&(jj>=0)&&(jj<w))
                                {
                                    if(ptr[jj])
                                        acum += (float)ptr[jj] *(float)((float)kernel_body_nm/(float)255.0f);
                                }
                            }
                            jj++;
                        }
                        ii++;
                    }
                    c_var = (int)(255.0*acum)/div;
                }
                else
                {
                    for(int n=0; n<kernel->h; n++)
                    {
                        int jj = j - dx;
                        const int * fast kernel_body_n = kernel->body + n*kernel->w;
                        for(int m=0; m<kernel->h; m++)
                        {
                            const int kernel_body_nm = kernel_body_n[m];
                            if(kernel_body_nm)
                            {
                                if((ii>=0)&&(ii<h)&&(jj>=0)&&(jj<w))
                                {
                                    if(ptr[jj])
                                        c_var += ptr[jj] * kernel_body_nm;
                                }
                            }
                            jj++;
                        }
                        ii++;
                    }
                    c_var/=div;
                }
                local->body[i*w+j] = c_var;
            }
        }
        return local;

    }

    void normalize()
    {
        int max_val = get(this->max());
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
                ptr[j] = 255*(float)(ptr[j])/(float)max_val;
        }
    }

    void threshold(int limit,bool binary=true)
    {
        for(int i=0; i<h; i++)
        {
            int * fast ptr = body+i*w;
            for(int j=0; j<w; j++)
            {
                if(ptr[j]<limit)
                    ptr[j]=0;
                else if(binary)
                {
                    ptr[j]=255;
                }
            }
        }
    }

    matrix * sobol_filterX()
    {
        matrix * kernel_mat = new matrix(3,3);
        int *fast kernel = kernel_mat->body;
        int kw = kernel_mat->w;
        kernel[0] = -1;
        kernel[1] =  0;
        kernel[2] =  1;
        kernel[kw] = -2;
        kernel[kw+1] =  0;
        kernel[kw+2] =  2;
        kernel[2*kw] = -1;
        kernel[2*kw+1] =  0;
        kernel[2*kw+2] =  1;

        matrix * local = conv(kernel_mat,false);

        delete kernel_mat;

        return local;
    }

    matrix * sobol_filterY()
    {
        matrix * kernel_mat = new matrix(3,3);
        int * fast kernel = kernel_mat->body;
        int kw = kernel_mat->w;
        kernel[0] = -1;
        kernel[1] = -2;
        kernel[2] = -1;
        kernel[kw] =  0;
        kernel[kw+1] =  0;
        kernel[kw+2] =  0;
        kernel[2*kw] =  1;
        kernel[2*kw+1] =  2;
        kernel[2*kw+2] =  1;

        matrix * local = conv(kernel_mat,false);

        delete kernel_mat;

        return local;
    }

    void sobol_filter()
    {
        matrix * gx = sobol_filterX();
        matrix * gy = sobol_filterX();
        matrix * g = new matrix(gx->w,gx->h);


        for(int i=0; i<gx->h; i++)
        {
            int * fast ptr = g->body+i*w;
            int * fast ptrx = gx->body+i*w;
            int * fast ptry = gy->body+i*w;
            for(int j=0; j<gx->w; j++)
            {
                const int gx_val = ptrx[j];
                const int gy_val = ptry[j];
                ptr[j] = sqrt((double)(gx_val*gx_val+gy_val*gy_val));
            }
        }


        delete gx;
        delete gy;
        set(*g);
        delete g;
    }

    void gaussian_filter(float sigma)
    {
        const int n = 2 * (int)(2 * sigma) + 3;
        const float mean = (float)floor(n / 2.0);
        matrix * kernel = new matrix(n,n);

        for (int i = 0; i < n; i++)
        {
            int * fast ptr = kernel->body+i*w;
            for (int j = 0; j < n; j++)
            {
                ptr[j] =(int)(255.0 * ( exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                                    pow((j - mean) / sigma, 2.0)))
                                        / (2 * M_PI * sigma * sigma)));
            }
        }

        matrix * local = conv(kernel,true);

        delete kernel;
        set(*local);
        delete local;
    }


};


struct SliceMatrix : public matrix
{
    size_t start_x;
    size_t start_y;
    SliceMatrix(matrix *m,size_t start_x,size_t start_y,size_t stop_x,size_t stop_y):matrix()
    {
        if(start_x>stop_x)
            std::swap(start_x,stop_x);
        if(start_y>stop_y)
            std::swap(start_y,stop_y);
        this->w = stop_x - start_x;
        this->h = stop_y - start_y;
        this->start_x = start_x;
        this->start_y = start_y;
        this->body = new int[w*h];
        for(int i=0; i<h; i++)
        {
            int * fast dest = body+i*w;
            int * fast src = m->body+(start_y+i)*w;
            for(int j=0; j<w; j++)
                dest[j] = src[start_x+j];
        }
    }

    void expand(size_t w,size_t h)
    {
        int *local = new int[w*h];
        for(size_t i=0; i<h; i++)
        {
            int * fast ptr = local+i*w;
            for(size_t j=0; j<w; j++)
                ptr[j] =0;
        }
        for(int ii=this->start_y,i=0; i<this->h; ii++,i++)
        {
            int * fast src = body+i*w;
            int * fast dest = local+ii*w;
            for(int jj=this->start_x,j=0; j<this->w; jj++,j++)
            {
                dest[jj] = src[j];
            }
        }
        delete[] body;
        body = local;
        this->w = w;
        this->h = h;
    }

};

#define BGR_BUFFER 1
#define RGB_BUFFER 0

class frame
{
private:
    int w;
    int h;
    matrix** rgb;
public:
    frame(int w,int h)
    {
        rgb = new matrix*[3];
        this->w=w;
        this->h=h;
        for(int i=0; i<3; i++)
            rgb[i]= new matrix(w,h);
    }
    frame(frame &f)
    {
        w = f.getWidth();
        h = f.getHeight();
        rgb = new matrix*[3];
        rgb[0] = new matrix(*f.getRedBuffer());
        rgb[1] = new matrix(*f.getGreenBuffer());
        rgb[2] = new matrix(*f.getBlueBuffer());


    }

    frame(matrix &m)
    {
        rgb = new matrix*[3];
        w = m.w;
        h = m.h;
        for(int i=0; i<3; i++)
            rgb[i] = new matrix(m);
    }
    frame(int *BUFFER,int w,int h,int RGB_OR_BGR=0):frame(w,h)
    {
        unsigned char * pixel;
        int r_chanel,b_chanel;
        if(RGB_OR_BGR)
        {
            r_chanel = 2;
            b_chanel = 0;
        }
        else
        {
            r_chanel = 0;
            b_chanel = 2;
        }

        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                pixel=(unsigned char*)&BUFFER[i*w+j];
                ptr_1[j]=pixel[r_chanel];
                ptr_2[j]=pixel[1];
                ptr_3[j]=pixel[b_chanel];
            }
        }
    }

    frame(Color *colors,int w,int h):frame(w,h)
    {
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                unsigned char * fast pixel=(unsigned char*)&colors[i*w+j];
                ptr_1[j]=pixel[0];
                ptr_2[j]=pixel[1];
                ptr_3[j]=pixel[2];
            }
        }
    }

        frame(const RGBImage & image):frame(image.width,image.height)
        {

            for(int i=0;i<h;i++)
            {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
                for(int j=0;j<w;j++)
                {
                    const int index =3*(i*w+j);
                    ptr_1[j]=image.data[index];
                    ptr_2[j]=image.data[index+1];
                    ptr_3[j]=image.data[index+2];
                }
            }

        }

    Color * toColorBuffer()
    {
        Color * pixel_buffer = new Color[h*w];
        for(int i =0 ; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                const int index = i*w+j;
                pixel_buffer[index].r = ptr_1[j];
                pixel_buffer[index].g = ptr_2[j];
                pixel_buffer[index].b = ptr_3[j];
                pixel_buffer[index].a = 0xFF;
            }

        }
        return pixel_buffer;
    }

    frame * subFrame(Point start,Point stop)
    {
        frame * local = new frame(stop.x-start.x,stop.y-start.y);

        int l_w = local->getWidth(), l_h = local->getHeight();

        for(int i=0; i<l_h; i++)
        {
            int * fast src_1 = rgb[0]->body+(start.y+i)*w;
            int * fast src_2 = rgb[1]->body+(start.y+i)*w;
            int * fast src_3 = rgb[2]->body+(start.y+i)*w;
            int * fast dest_1 = local->rgb[0]->body+i*l_w;
            int * fast dest_2 = local->rgb[1]->body+i*l_w;
            int * fast dest_3 = local->rgb[2]->body+i*l_w;

            for(int j=0; j<l_w; j++)
            {
                dest_1[j] = src_1[start.x+j];
                dest_2[j] = src_2[start.x+j];
                dest_3[j] = src_3[start.x+j];
            }
        }
        return local;
    }

    ~frame()
    {
        for(int i=0; i<3; i++)
            delete rgb[i];
        delete[] rgb;
    }

    int getWidth()
    {
        return w;
    }

    int getHeight()
    {
        return h;
    }

    const matrix* getRedBuffer()
    {
        return rgb[0];
    }
    const matrix* getGreenBuffer()
    {
        return rgb[1];
    }
    const matrix* getBlueBuffer()
    {
        return rgb[2];
    }

    frame* clone()
    {
        return new frame(*this);
    }

    matrix* toGreyScaleBuffer()
    {
        matrix *m = new matrix(w,h);
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            int * fast ptr = m->body+i*w;
            for(int j=0; j<w; j++)
            {
                ptr[j] = 0.4*ptr_1[j]+
                                0.3*ptr_2[j]+
                                0.3*ptr_3[j];
            }
        }
        return m;
    }

    Color getPixel(int x,int y)
    {
        Color c;
        c.a=0xff;
        c.r = rgb[0]->body[y*w+x];
        c.g = rgb[1]->body[y*w+x];
        c.b = rgb[2]->body[y*w+x];
        return c;
    }

    void setPixel(int x,int y,Color c)
    {
        rgb[0]->body[y*w+x] = c.r;
        rgb[1]->body[y*w+x] = c.g;
        rgb[2]->body[y*w+x] = c.b;
    }

    matrix** getBuffers()
    {
        return rgb;
    }

    void setPixelBuffer(uint8_t *pixel_buff,bool RGB = true)
    {

        int red = (RGB)?0:2;
        int blue = (RGB)?2:0;
        int k = 0 ;
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[red]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[blue]->body+i*w;
            for(int j=0; j<w; j++)
            {
                pixel_buff[k] = ptr_1[j];
                pixel_buff[k+1] = ptr_2[j];
                pixel_buff[k+2] = ptr_3[j];
                k+=3;
            }
        }
    }

    void reduce_white(float error,bool toMax=false)
    {
        float loss = 0.0f;
        if(error>1.1f)
            error/=100.0f;
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];

                loss = sqrt(pow((*r)/(float)(*g)-1,2)+pow((*g)/(float)(*b)-1,2));
                if(loss>error)
                {
                    *r = 0;
                    *g = 0;
                    *b = 0;
                }
                else if(toMax)
                {
                    *r = 255;
                    *g = 255;
                    *b = 255;
                }

            }
        }

    }

    void reduce_HSV(float angle,float saturation,float value,float error,bool toMax = false)
    {
        if(error>1.1f)
            error/=100.0f;
        if(saturation>1.1f)
            saturation/=100.0f;
        if(value>1.1f)
            value/=100.0f;

        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];
                Color c = {(unsigned char)*r,(unsigned char)*g,(unsigned char)*b,(unsigned char)0xFF};
                Vector3 a = ColorToHSV(c);
                //printf("%f \n",a.x);
                a.x-=angle;
                a.y-=saturation;
                a.z-=value;
                //loss = 0.3*sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
                if((fabs(a.x)>error)||(fabs(a.y)>error)||(fabs(a.z)>error))
                {
                    *r = 0;
                    *g = 0;
                    *b = 0;
                }
                else if(toMax)
                {
                    *r = 255;
                    *g = 255;
                    *b = 255;
                }
            }
        }
    }
    void reduce_red(float error,bool toMax = false)
    {

        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];
                Color c = {(unsigned char)*r,(unsigned char)*g,(unsigned char)*b,0xFF};
                Vector3 a = ColorToHSV(c);

                if((abs(355-a.x)>error)||((a.y<0.7)&&(a.z<0.7)))
                {
                    *r = 0;

                }
                else if(toMax)
                {
                    *r = 255;
                }

                *g = *r;
                *b = *r;
            }
        }
    }

    void reduce_blue()
    {

        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];
                int db = (*r+*g)/2, dg = (*r+*b)/2,dr =(*b+*g)/2;
                if(*r<dr)
                    *r=0;
                else *r -= dr;
                if(*g<dg)
                    *g =0;
                else *g -= dg;
                if(*b<db)
                    *b =0;
                else *b -= db;
            }
        }
    }

    void reduce_space(float saturation,float value,float error,bool toMax = false)
    {
        float loss =0.0f;
        if(error>1.1f)
            error/=100.0f;
        if(saturation>1.1f)
            saturation/=100.0f;
        if(value>1.1f)
            value/=100.0f;
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];
                Color c = {(unsigned char)*r,(unsigned char)*g,(unsigned char)*b,0xFF};
                Vector3 a = ColorToHSV(c);
                a.y-=saturation;
                a.z-=value;
                loss = sqrt(a.y*a.y + a.z*a.z);
                if(loss>error)
                {
                    *r = 0;
                    *g = 0;
                    *b = 0;
                }
                else if(toMax)
                {
                    *r = 255;
                    *g = 255;
                    *b = 255;
                }
            }
        }
    }

    void reduce_spectre()
    {
        rgb[0]->threshold(150,false);
        rgb[1]->threshold(150,false);
        rgb[2]->threshold(150,false);
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];
                Color c = {(unsigned char)*r,(unsigned char)*g,(unsigned char)*b,0xFF};
                Vector3 a = ColorToHSV(c);
                if(a.x>5){
                //a.y=1.0f;
                a.z=1.0f;
                c = ColorFromHSV(a);
                *r =c.r;
                *g =c.g;
                *b =c.b;
                }
            }
        }
    }

    void balance_space(float saturation,float value)
    {
        if(saturation>1.1f)
            saturation/=100.0f;
        if(value>1.1f)
            value/=100.0f;
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                int * r = &ptr_1[j];
                int * g = &ptr_2[j];
                int * b = &ptr_3[j];
                Color c = {(unsigned char)*r,(unsigned char)*g,(unsigned char)*b,0xFF};
                Vector3 a = ColorToHSV(c);
                a.y=saturation;
                a.z=value;
                c = ColorFromHSV(a);
                *r = c.r;
                *g = c.g;
                *b = c.b;

            }
        }
    }

    void normalize()
    {

        for(int i=0; i<3; i++)
            rgb[i]->normalize();
    }

    void setColorBuffer(Color * pixel_buff)
    {
        for(int i=0; i<h; i++)
        {
            int * fast ptr_1 = rgb[0]->body+i*w;
            int * fast ptr_2 = rgb[1]->body+i*w;
            int * fast ptr_3 = rgb[2]->body+i*w;
            for(int j=0; j<w; j++)
            {
                Color c;
                c.r = ptr_1[j];
                c.g = ptr_2[j];
                c.b = ptr_3[j];
                c.a = 0xFF;
                pixel_buff[i*w+j]=c;
            }
        }
    }

    void setSize(int w,int h)
    {
        this->w = w;
        this->h = h;
        for(int i=0; i<3; i++)
            rgb[i]->setSize(w,h);
    }

    void scale(float ration)
    {
        this->w *= ration;
        this->h *= ration;
        for(int i=0; i<3; i++)
            rgb[i]->setSize(this->w,this->h);
    }


};



