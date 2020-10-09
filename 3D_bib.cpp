
#include "3D_bib.h"
#include <iostream>

using namespace std;

Operaciones3D::Operaciones3D()
{
        //Variables para operaciones trigonometricas
        pi = 3.14159265359;

}


//recordar que (pi/180 = r/g) donde "r" son radianes y "g" grados
//se aplica la formula r
float Operaciones3D::RadToDeg(float r)
{
      return ((180*r)/pi);
}

float Operaciones3D::DegToRad(float g)
{
      return ((g*pi)/180);
}

void Operaciones3D::LoadIdentity(float M[][4])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      if(i==j)
         M[i][j]=1;
      else
         M[i][j]=0;
 }

void Operaciones3D::translate(float x, float y, float z)
{
  LoadIdentity(T);
  T[0][3]=x;
  T[1][3]=y;
  T[2][3]=z;
 }

void Operaciones3D::rotateX(float deg)
{
  LoadIdentity(R);
  R[1][1] = cos(deg);
  R[1][2] = -1*sin(deg);
  R[2][1] = sin(deg);
  R[2][2] = cos(deg);
 }

void Operaciones3D::rotateY(float deg)
{
  LoadIdentity(R);
  R[0][0] = cos(deg);
  R[0][2] = sin(deg);
  R[2][0] = -1*sin(deg);
  R[2][2] = cos(deg);
 }

void Operaciones3D::rotateZ(float deg)
{
  LoadIdentity(R);
  R[0][0] = cos(deg);
  R[0][1] = -1*sin(deg);
  R[1][0] = sin(deg);
  R[1][1] = cos(deg);
 }

void Operaciones3D::MultM(float M1[][4], float M2[][4], float Res[][4])
{
  float tmp[4][4];
  int i,j,k;
  for(i=0; i<4;i++)
     for(j=0;j<4;j++){
        tmp[i][j]=0;
        for(k=0;k<4;k++)
           tmp[i][j]+=M1[i][k]*M2[k][j];
     }
  for(i=0;i<4;i++)
     for(j=0;j<4;j++)
        Res[i][j] = tmp[i][j];
}

//multiplica la matriz m por el punto p y regresa el resultado en el punto p
void Operaciones3D::MatPoint(float m[][4], float p[3])
{
  float tmp[4];
  int i,j;
  for(i=0; i<3; i++)
    { tmp[i] = p[i];
      p[i] = 0;
    }
  tmp[3]=1;
  for(i=0;i<3;i++)
    for(j=0;j<4;j++)
        p[i] += m[i][j]*tmp[j];
}

//multiplica la matriz m por cada punto del objeto definido por la matriz p de size x 3
void Operaciones3D::MatObject(float m[][4], int size, float p[][3])
{
     int i;
     for(i=0; i<size; i++)
       MatPoint(m,p[i]);
}

//rotacion paralela a uno de los ejes
//theta: angulo de rotacion;
//distA,distB: vector (distA,distB) que separa al eje de rotacion del objeto
//con respecto a uno de los ejes del sistema carteciano. Si el eje es:
//X: (distA,distB) es el vector (0,distA,distB)
//Y: (distA,distB) es el vector (distA,0,distB)
//Z: (distA,distB) es el vector (distA,distB,0)
void Operaciones3D::RotacionParalela(char eje, float theta, float distA, float distB)
{
     switch(eje){
        case 'X':
             //se actualiza la matriz de traslacion para mover el objeto en el espacio
             translate(0,-distA,-distB);
             //se actualiza la matriz de rotacion con el angulo especificado
             rotateX(DegToRad(theta));
             //se multiplica la matriz de rotacion por la de traslacion actual
             //el resultado queda en la matriz A
             MultM(R,T,A);
             //se agrega la matriz de traslacion inversa a A
             translate(0,distA,distB);
             //se multiplica la matriz de traslacion por la matriz A y se deja el resultado en A
             MultM(T,A,A);
             break;
        case 'Y':
            //se actualiza la matriz de traslacion para mover el objeto en el espacio
             translate(-distA,0,-distB);
             //se actualiza la matriz de rotacion con el angulo especificado
             rotateY(DegToRad(theta));
             //se multiplica la matriz de rotacion por la de traslacion actual
             //el resultado queda en la matriz A
             MultM(R,T,A);
             //se agrega la matriz de traslacion inversa a A
             translate(distA,0,distB);
             //se multiplica la matriz de traslacion por la matriz A y se deja el resultado en A
             MultM(T,A,A);
             break;
        case 'Z':
            //se actualiza la matriz de traslacion para mover el objeto en el espacio
             translate(-distA,-distB,0);
             //se actualiza la matriz de rotacion con el angulo especificado
             rotateZ(DegToRad(theta));
             //se multiplica la matriz de rotacion por la de traslacion actual
             //el resultado queda en la matriz A
             MultM(R,T,A);
             //se agrega la matriz de traslacion inversa a A
             translate(distA,distB,0);
             //se multiplica la matriz de traslacion por la matriz A y se deja el resultado en A
             MultM(T,A,A);
             break;
     }
}

float getNorma(float a, float b, float c)
{
    return sqrtf(pow(a,2)+pow(b,2)+pow(c,2));
}

float getNorma(float b, float c)
{
    return sqrtf(pow(b,2)+pow(c,2));
}

void Operaciones3D::rotateAlpha(float b, float c, float d, bool inv)
{
    LoadIdentity(R);
    R[1][1] = c/d;
    R[2][1] = inv ? -b/d : b/d;
    R[1][2] = inv ? b/d : -b/d;
    R[2][2] = c/d;
}

void Operaciones3D::rotateBeta(float a, float d, bool inv)
{
    LoadIdentity(R);
    R[0][0] = d;
    R[0][2] = inv ? -a : a;
    R[2][0] = inv ? a : -a;
    R[2][2] = d;
}

void Operaciones3D::RotacionLibre(float theta, float p1[3], float p2[3])
{
    float x = p2[0]-p1[0];
    float y = p2[1]-p1[1];
    float z = p2[2]-p1[2];
    float v = getNorma(x,y,z);
    float a = x/v;
    float b = y/v;
    float c = z/v;
    float d = getNorma(b,c);
    translate(-p1[0],-p1[1],-p1[2]);
    rotateAlpha(b,c,d,false);
    MultM(R,T,A);

    rotateBeta(a,d,false);
    MultM(R,A,A);

    rotateZ(DegToRad(theta));
    MultM(R,A,A);

    rotateBeta(a,d,true);
    MultM(R,A,A);

    rotateAlpha(b,c,d,true);
    MultM(R,A,A);

    translate(p1[0],p1[1],p1[2]);
    MultM(T,A,A);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
