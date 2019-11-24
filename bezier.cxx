#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <cmath>

using std::cerr;
using std::endl;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

double dot_product(double* A, double* B)
{
    return (A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2]);
}

double* cross_product(double* A, double* B)
{
    double* out = new double[3];
    out[0] = A[1] * B[2] - A[2] * B[1];
    out[1] = A[2] * B[0] - A[0] * B[2];
    out[2] = A[0] * B[1] - A[1] * B[0];

    return out;
}

double* normalize(double x, double y, double z)
{
    double* out = new double[3];
    double length = sqrt( (x*x + y*y + z*z));
    out[0] = x / length;
    out[1] = y / length;
    out[2] = z / length;
    return out;
}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    double*         u;
    double*         v;
    double*         w;

    // use your own project 1 code for those three functions
    void            positionSetup(void);
    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(int n, int m);
};

void Camera::positionSetup(void)
{
    double* new_up = new double[3];
    u = new double[3];
    v = new double[3];
    w = new double[3];

    for(int i = 0; i < 3; i++)
    {
        w[i] = position[i] - focus[i];
        new_up[i] = up[i];
    }

    u = cross_product(new_up, w);
    v = cross_product(w, u);

    u = normalize(u[0], u[1], u[2]);
    v = normalize(v[0], v[1], v[2]);
    w = normalize(w[0], w[1], w[2]);

    delete new_up;

}

Matrix Camera::ViewTransform(void)
{
    Matrix Transform_view;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Transform_view.A[i][j] = 0;
        }
    }

    Transform_view.A[0][0] = 1 / (tan(angle / 2));

    Transform_view.A[1][1] = 1 / (tan(angle / 2));

    Transform_view.A[2][2] = (far + near) / (far - near);
    Transform_view.A[2][3] = -1;

    Transform_view.A[3][2] = (2 * far * near) / (far - near);

    return Transform_view;
}

Matrix Camera::CameraTransform(void)
{
    Matrix Transform_camera;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Transform_camera.A[i][j] = 0;
        }
    }

    double* t = new double[3];
    t[0] = 0 - position[0];
    t[1] = 0 - position[1];
    t[2] = 0 - position[2];

    Transform_camera.A[0][0] = u[0];  //u.x
    Transform_camera.A[0][1] = v[0];  //v.x
    Transform_camera.A[0][2] = w[0];  //w.x

    Transform_camera.A[1][0] = u[1];  //u.y
    Transform_camera.A[1][1] = v[1];  //v.y
    Transform_camera.A[1][2] = w[1];  //w.y

    Transform_camera.A[2][0] = u[2];  //u.z
    Transform_camera.A[2][1] = v[2];  //v.z
    Transform_camera.A[2][2] = w[2];  //w.z

    Transform_camera.A[3][0] = dot_product(u, t);  //u dot t
    Transform_camera.A[3][1] = dot_product(v, t);  //v dot t
    Transform_camera.A[3][2] = dot_product(w, t);  //w dot t
    Transform_camera.A[3][3] = 1;

    delete t;
    return Transform_camera;
}

Matrix Camera::DeviceTransform(int n, int m)
{
    Matrix Transform_Device;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Transform_Device.A[i][j] = 0;
        }
    }

    double x = n / 2;
    double y = m / 2;
    Transform_Device.A[0][0] = x;

    Transform_Device.A[1][1] = y;

    Transform_Device.A[2][2] = 1;

    Transform_Device.A[3][0] = x;
    Transform_Device.A[3][1] = y;
    Transform_Device.A[3][3] = 1;

    return Transform_Device;
}

// note that this is different from project 1 starter code
Camera
GetCamera(int frame, int nframes)
{
    double t = (double)frame/(double)nframes;
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 80 * sin(2*M_PI*t);
    c.position[1] = 30;
    c.position[2] = 80 * cos(2*M_PI*t);
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

class Screen
{
    public:
        unsigned char   *buffer;
        double          *zbuffer;
        int width, height;

        void colorPixel(int r, int c, unsigned char color[3], double depth);

  // would some methods for accessing and setting pixels be helpful?
};

void Screen::colorPixel(int r, int c, unsigned char color[3], double depth){
    if (r < 0 || r >= height || c < 0 || c >= width) return;

    int index = r * width + c;

    if (depth > zbuffer[index])
    {
        buffer[3*index+0] = color[0];
        buffer[3*index+1] = color[1];
        buffer[3*index+2] = color[2];
        zbuffer[index] = depth;
    }
}

class Line{
public:
    double X[2];
    double Y[2];
    double Z[2];
    unsigned char   color[3];

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }
    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Draw(Screen * screen){
        int x0, y0, x1, y1;
        if(Y[0] < Y[1]){
            x0 = X[0];
            y0 = Y[0];
            x1 = X[1];
            y1 = Y[1];
        }else{
            x0 = X[1];
            y0 = Y[1];
            x1 = X[0];
            y1 = Y[0];
        }

        int e;
        int dx = x1 - x0;
        int dy = y1 - y0;
        int x = x0;
        int y = y0;
        double dz = Z[1] - Z[0];

        int stepx = 1;
        int stepy = 1;
        if (dx<0){
            dx = - dx;
            stepx = -1;
        }
        if (dy<0){
            dy = - dy;
            stepy = -1;
        }

        int i;
        // first octant
        if(dy<dx){
            e = 2 * dy - dx;
            for(i=0; i<dx; i++){
                screen->colorPixel(y, x, color, i*dz/dx+Z[0]);
                while(e > 0){
                    y += stepy;
                    e = e - 2 * dx;
                }
                x += stepx;
                e = e + 2 * dy;
            }
        // second octant
        }else{
            e = 2 * dx - dy;
            for(i=0; i<dy; i++){
                screen->colorPixel(y, x, color, i*dz/dy+Z[0]);
                while(e > 0){
                    x += stepx;
                    e = e - 2 * dy;
                }
                y += stepy;
                e = e + 2 * dx;
            }
        }
    }

    void Print(){
        std::cout << "v1: " << X[0] << ", " << Y[0] << ", " << Z[0] << '\n';
        std::cout << "v2: " << X[1] << ", " << Y[1] << ", " << Z[1] << '\n';
    }
};

void Casteljau(double p[4], double q[4], double r[4])
{
    q[0] = p[0];
    q[1] = 0.5 * (p[0] + p[1]);
    q[2] = 0.5 * q[1] + 0.25 * (p[1] + p[2]);
    r[2] = 0.5 * (p[2] + p[3]);
    r[1] = 0.5 * r[2] + 0.25 * (p[1] + p[2]);
    q[3] = 0.5 * (q[2] + r[1]);
    r[3] = p[3];
    r[0] = q[3];
}

void Bezier_Divide(double pX[4], double pY[4], double pZ[4], double qX[4], double qY[4], double qZ[4], double rX[4], double rY[4], double rZ[4]){
    // your code goes here
    Casteljau(pX, qX, rX);
    Casteljau(pY, qY, rY);
    Casteljau(pZ, qZ, rZ);
}

class BezierCurve{
public:
    double X[4];
    double Y[4];
    double Z[4];
    unsigned char color[3];

    double* ptIn;
    double* ptOut;

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }

    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Print(){
        printf("p0: %lf, %lf, %lf\n", X[0], Y[0], Z[0]);
        printf("p1: %lf, %lf, %lf\n", X[1], Y[1], Z[1]);
        printf("p2: %lf, %lf, %lf\n", X[2], Y[2], Z[2]);
        printf("p3: %lf, %lf, %lf\n", X[3], Y[3], Z[3]);
    }

    void RotateY(BezierCurve * dst, double angleDegree){
        int i;
        for(i=0; i<4; i++){
            dst->X[i] = X[i] * cos(2*M_PI*angleDegree/360)
                        + Z[i] * sin(2*M_PI*angleDegree/360);
            dst->Y[i] = Y[i];
            dst->Z[i] = Z[i] * cos(2*M_PI*angleDegree/360)
                        - X[i] * sin(2*M_PI*angleDegree/360);
        }
    }

    void curve_transform(Matrix m) {
        for (int i = 0; i < 4; i++)
        {
            ptIn = new double[4];
            ptOut = new double[4];

            ptIn[0] = X[i];
            ptIn[1] = Y[i];
            ptIn[2] = Z[i];
            ptIn[3] = 1;

            m.TransformPoint(ptIn, ptOut);

            ptOut[0] = ptOut[0] / ptOut[3];
            ptOut[1] = ptOut[1] / ptOut[3];
            ptOut[2] = ptOut[2] / ptOut[3];

            X[i] = ptOut[0];
            Y[i] = ptOut[1];
            Z[i] = ptOut[2];

        }
    }
    
    void Draw(Screen * screen){
        double limit = 500;
        double qx[4], qy[4], qz[4], rx[4], ry[4], rz[4];
        double current, max;
        current = X[0] - X[1];
        current = max = current * current;
        for (int i = 0; i < 3; i++)
        {
            for (int j = i+1; j<4; j++){
                current = X[i] - X[i+1];
                current = current * current;
                if (current > max) {
                    max = current;
                }
                current = Y[i] - Y[i+1];
                current = current * current;
                if (current > max) {
                    max = current;
                }
            }
        }
        Bezier_Divide(X, Y, Z, qx, qy, qz, rx, ry, rz);
        if (max < limit) {
            for (int i = 0; i < 3; i++)
            {
                Line l1, l2;
                l1.SetPointByIndex(0, qx[i], qy[i], qz[i]);
                l1.SetPointByIndex(1, qx[i+1], qy[i+1], qz[i+1]);

                l2.SetPointByIndex(0, rx[i], ry[i], rz[i]);
                l2.SetPointByIndex(1, rx[i+1], ry[i+1], rz[i+1]);

                l1.SetColor(color[0], color[1], color[2]);
                l2.SetColor(color[0], color[1], color[2]);

                l1.Draw(screen);
                l2.Draw(screen);
            }
        } else {
            BezierCurve c1, c2;
            for (int i = 0; i < 4; i++)
            {
                c1.SetPointByIndex(i, qx[i], qy[i], qz[i]);
                c2.SetPointByIndex(i, rx[i], ry[i], rz[i]);
            }
            c1.SetColor(color[0], color[1], color[2]);
            c2.SetColor(color[0], color[1], color[2]);

            c1.Draw(screen);
            c2.Draw(screen);

        }
    }
};

class BezierSurface{
public:
    double X[16];
    double Y[16];
    double Z[16];
    unsigned char color[3];

    double qx[4], qy[4], qz[4], rx[4], ry[4], rz[4];
    double sx[16], sy[16], sz[16], srx[16], sry[16], srz[16];

    double* ptIn;
    double* ptOut;

    void SetPoints(double x[16], double y[16], double z[16]){
        int i;
        for(i=0; i<16; i++){
            X[i] = x[i];
            Y[i] = y[i];
            Z[i] = z[i];
        }
    }

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }

    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void surface_transform(Matrix m)
    {  
        for (int i = 0; i < 16; i++)
        {
            ptIn = new double[4];
            ptOut = new double[4];

            ptIn[0] = X[i];
            ptIn[1] = Y[i];
            ptIn[2] = Z[i];
            ptIn[3] = 1;

            m.TransformPoint(ptIn, ptOut);

            ptOut[0] = ptOut[0] / ptOut[3];
            ptOut[1] = ptOut[1] / ptOut[3];
            ptOut[2] = ptOut[2] / ptOut[3];

            X[i] = ptOut[0];
            Y[i] = ptOut[1];
            Z[i] = ptOut[2];
        }
    }

    void Print(){
        printf("v0: %lf, %lf, %lf\n", X[0], Y[0], Z[0]);
        printf("v1: %lf, %lf, %lf\n", X[1], Y[1], Z[1]);
        printf("v2: %lf, %lf, %lf\n", X[2], Y[2], Z[2]);

    }

    void Set_Curve(BezierCurve c, double x[4], double y[4], double z[4], Screen* s)
    {
        for (int i = 0; i < 4; i++)
        {
            c.SetPointByIndex(i, x[i], y[i], z[i]);
        }
        c.SetColor(color[0], color[1], color[2]);
        c.Draw(s);
    }

    void subsurface(BezierSurface bs, Screen* s, int divisions)
    {
        BezierSurface s1, s2;
        for (int r = 0; r < 4; r++)
        {
            double x[4], y[4], z[4];
            for (int c = 0; c < 4; c++)
            {
                x[c] = bs.X[4*c+r];
                y[c] = bs.Y[4*c+r];
                z[c] = bs.Z[4*c+r];
            }
            Bezier_Divide(x, y, z, qx, qy, qz, rx, ry, rz);
            for (int c = 0; c < 4; c++)
            {
                sx[4*r+c] = qx[c];
                sy[4*r+c] = qy[c];
                sz[4*r+c] = qz[c];
                srx[4*r+c] = rx[c];
                sry[4*r+c] = ry[c];
                srz[4*r+c] = rz[c];   
            }
        }
        s1.SetPoints(sx, sy, sz);
        s2.SetPoints(srx, sry, srz);

        s1.SetColor(color[0], color[1], color[2]);
        s2.SetColor(color[0], color[1], color[2]);

        s1.Draw(s, divisions);
        s2.Draw(s, divisions);
    }

    void Draw(Screen * screen, int divisions){
        if (divisions == 0)
        {
            for (int r = 0; r < 4; r++)
            {
                double qx[4], qy[4], qz[4], rx[4], ry[4], rz[4];
                double rowx[4], rowy[4], rowz[4];
                double colx[4], coly[4], colz[4];

                for (int c = 0; c < 4; c++)
                {
                    rowx[c] = X[4*r+c];
                    rowy[c] = Y[4*r+c];
                    rowz[c] = Z[4*r+c];
                    colx[c] = X[4*c+r];
                    coly[c] = Y[4*c+r];
                    colz[c] = Z[4*c+r];
                }

                BezierCurve c1, c2, c3, c4;
                Bezier_Divide(rowx, rowy, rowz, qx, qy, qz, rx, ry, rz);
                Set_Curve(c1, qx, qy, qz, screen);
                Set_Curve(c2, rx, ry, rz, screen);
                Bezier_Divide(colx, coly, colz, qx, qy, qz, rx, ry, rz);
                Set_Curve(c3, qx, qy, qz, screen);
                Set_Curve(c4, rx, ry, rz, screen);

            }
        }
        else 
        {
            divisions--;
            BezierSurface s1, s2;
            for (int r = 0; r < 4; r++)
            {
                double x[4], y[4], z[4];
                for (int c = 0; c < 4; c++)
                {
                    x[c] = X[4*r+c];
                    y[c] = Y[4*r+c];
                    z[c] = Z[4*r+c];
                }
                Bezier_Divide(x, y, z, qx, qy, qz, rx, ry, rz);
                for (int c = 0; c < 4; c++)
                {
                    sx[4*r+c] = qx[c];
                    sy[4*r+c] = qy[c];
                    sz[4*r+c] = qz[c];
                    srx[4*r+c] = rx[c];
                    sry[4*r+c] = ry[c];
                    srz[4*r+c] = rz[c];   
                }
            }
            s1.SetPoints(sx, sy, sz);
            s2.SetPoints(srx, sry, srz);
            subsurface(s1, screen, divisions);
            subsurface(s2, screen, divisions);
        }
        
    }
};

// returns an array of BezierCurves of size 6
BezierCurve * getLeaves(){
    BezierCurve *c = (BezierCurve*)malloc(sizeof(BezierCurve)*6);
    int i;
    for(i=0; i<6; i++)
        c[i].SetColor(30, 215, 97);

    c[0].SetPointByIndex(0, 0, -10, 0);
    c[0].SetPointByIndex(1, 14, -7.2, 0);
    c[0].SetPointByIndex(2, 9.8, -3, 2.8);
    c[0].SetPointByIndex(3, 14, 4, 14);

    c[1].SetPointByIndex(0, 0, -10, 0);
    c[1].SetPointByIndex(1, 0, -7.2, 14);
    c[1].SetPointByIndex(2, 2.8, -3, 9.8);
    c[1].SetPointByIndex(3, 14, 4, 14);

    c[0].RotateY(&(c[2]), 120);
    c[1].RotateY(&(c[3]), 120);
    c[0].RotateY(&(c[4]), 240);
    c[1].RotateY(&(c[5]), 240);

    return &(c[0]);
}

// returns an array of BezierSurfaces of size 2
BezierSurface * getSurfaces(){
    BezierSurface *s = (BezierSurface*)malloc(sizeof(BezierSurface)*2);

    s[0].SetPointByIndex(0, 0.0, 0.0, 0.0); // first row
	s[0].SetPointByIndex(1, 0.0, 3, 5);
	s[0].SetPointByIndex(2, 0.0, 7.5, 10);
	s[0].SetPointByIndex(3, 0.0, 1.5, 15.0);
	s[0].SetPointByIndex(4, 5, 12, 0.0); // second row
	s[0].SetPointByIndex(5, 5, -1.5, 5);
	s[0].SetPointByIndex(6, 5, 0.0, 10);
	s[0].SetPointByIndex(7, 5, 1.5, 15.0);
	s[0].SetPointByIndex(8, 10, 4.5, 0.0); // third row
	s[0].SetPointByIndex(9, 10, 12, 5);
	s[0].SetPointByIndex(10, 10, 13.5, 10);
	s[0].SetPointByIndex(11, 10, 7.5, 15.0);
	s[0].SetPointByIndex(12, 15.0, 6, 0.0); // fourth row
	s[0].SetPointByIndex(13, 15.0, 3, 5);
	s[0].SetPointByIndex(14, 15.0, 7.5, 10);
	s[0].SetPointByIndex(15, 15.0, 15.0, 15.0);
    s[0].SetColor(51, 133, 229);

    s[1].SetPointByIndex(0, 0.0, -3, 0.0); // first row
	s[1].SetPointByIndex(1, 0.0, -3, 5);
	s[1].SetPointByIndex(2, 0.0, -3, 10);
	s[1].SetPointByIndex(3, 0.0, -3, 15);
	s[1].SetPointByIndex(4, 5, -3, 0.0); // second row
	s[1].SetPointByIndex(5, 5, -3, 5);
	s[1].SetPointByIndex(6, 5, -3, 10);
	s[1].SetPointByIndex(7, 5, -3, 15);
	s[1].SetPointByIndex(8, 10, -3, 0.0); // third row
	s[1].SetPointByIndex(9, 10, -3, 5);
	s[1].SetPointByIndex(10, 10, -3, 10);
	s[1].SetPointByIndex(11, 10, -3, 15);
	s[1].SetPointByIndex(12, 15, -3, 0.0); // fourth row
	s[1].SetPointByIndex(13, 15, -3, 5);
	s[1].SetPointByIndex(14, 15, -3, 10);
	s[1].SetPointByIndex(15, 15, -3, 15);
    s[1].SetColor(31, 179, 83);

    return &(s[0]);
}

int main()
{
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer =
        (unsigned char *) image->GetScalarPointer(0,0,0);

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;
    screen.zbuffer = (double*)malloc(sizeof(double) * screen.width * screen.height);

    // uncooment the following two lines to get the curves and surfaces
    //BezierCurve *c = getLeaves();
    //BezierSurface * s = getSurfaces();
    // Camera rotates around y-axis
    // It comes back to the original position after 100 iterations
    int i;
    for(i=0; i<100; i++){
        int j;
        int npixels = screen.width * screen.height;
        for (int j = 0 ; j < npixels*3 ; j++)
            screen.buffer[j] = 0;
        for(int j = 0; j < npixels; j++)
            screen.zbuffer[j] = -1;

        Camera camera = GetCamera(i*10, 1000);
        camera.positionSetup();
        Matrix ct = camera.CameraTransform();
        Matrix vt = camera.ViewTransform();
        Matrix intermediate = Matrix::ComposeMatrices(ct, vt);
        Matrix dt = camera.DeviceTransform(screen.width, screen.height);
        Matrix total = Matrix::ComposeMatrices(intermediate, dt);

        BezierCurve *c = getLeaves();
        BezierSurface *s = getSurfaces();

        // for (int t = 0; t < 6; t++)
        // {
        //     c[t].curve_transform(total);
        //     c[t].Draw(&screen);
        // }  

        for (int t = 0; t < 2; t++)
        {
            s[t].surface_transform(total);
            s[t].Draw(&screen, 5);
        } 
        // draw your curves and surfaces here

        char name[20];
        // make sure you have a directory named images
        // so that writing images won't cause you errors
        sprintf(name, "./images/frame%02d", i);
        WriteImage(image, name);
    }
}