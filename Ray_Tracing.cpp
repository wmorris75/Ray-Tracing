/* SimpleRay main.c */
/* raytrace spheres and planes */

#include<stdio.h>
#include <GL/glut.h>              
#include <math.h>

// object type codes
#define PLANE 1        
#define SPHERE 2 
      
// screen size to raytrace
#define SCRNWIDTH 600       
#define SCRNHEIGHT 600       

#define NUMOBJECTS 6    
#define MAXLEVEL 4    
#define TOL 0.001      

typedef struct
{
  double x, y, z;
} VECTOR;

typedef struct
{       
  double r, g, b;       
} COLOR;               
                       
typedef struct
{       
  int Shape;  
  COLOR kd;
  double ks;    
  int Nzero;                          
  double A, B, C, D;    //sphere position,radius OR plane coeffiecients
} OBJECT;

typedef struct
{       
  VECTOR position;          
  COLOR Color;         
} LIGHT;

COLOR Iambient={0.2,0.2,0.2};
LIGHT Light = {{25.0, 170.0, 70.0},        // x, y, z
               {0.525, 0.525,0.525}};      // r, g, b

// This array contains the objects to display
OBJECT Objects[NUMOBJECTS] =
  {{SPHERE,                 //large shiny
    {0.15, 0.15, 0.15},
    0.9,
	20,
    2.0, 1.75, 3.0, 1.0},              
   {SPHERE,                 //little purple
    {0.7, 0.0, 0.9},
    0.1, 
    20,
    3.5, 0.75, 2.25, 0.25},              
   {PLANE,                  // grey floor
    {0.42, 0.42, 0.42},
    0.3, 
    4,
    0.0, 1.0, 0.0, 0.0},               
   {PLANE,                  // right
    {0.0, 0.0, 1.0},
    0.1,
    4,
    1.0, 0.0, 0.0, 0.0},               
   {SPHERE,                 // large dull red
    {1.0, 0.0, 0.0},
    0.1,
    1, 
    3.5, 1.0, 7.0, 0.75},               
   {PLANE,                  // left
    {0.2, 1.0, 0.0},
    0.1, 
    4, 
    0.0, 0.0, 1.0, 0.0}
  };              

COLOR background = {0, 0, 0};         

//VECTOR From, At, Up, A1, A2, A3;
//double DVal, VuAngle;
//double ImWd = SCRNWIDTH, ImHt = SCRNHEIGHT;

COLOR TraceRay(int Level, VECTOR Base, VECTOR Dir);

// some useful vector functions
VECTOR Subtract(VECTOR V1, VECTOR V2)
{
  VECTOR D;

  D.x = V1.x - V2.x;
  D.y = V1.y - V2.y;
  D.z = V1.z - V2.z;
  return D;
}

VECTOR Cross(VECTOR V1, VECTOR V2)
{
  VECTOR C;

  C.x = V1.y * V2.z - V2.y * V1.z;
  C.y = V1.z * V2.x - V2.z * V1.x;
  C.z = V1.x * V2.y - V2.x * V1.y;
  return C;
}

VECTOR Normalize(VECTOR V)
{
  VECTOR T;
  double D = sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
  
  T.x = V.x / D;
  T.y = V.y / D;
  T.z = V.z / D;
  return T;
}

double Dot(VECTOR V1, VECTOR V2)
{
  return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
}

// intersection functions
double IntersectPlane(OBJECT Obj, VECTOR Base, VECTOR Dir)
{
  double denom;

  denom = Obj.A * Dir.x + Obj.B * Dir.y + Obj.C * Dir.z;
  if (denom != 0)
    return -(Obj.A * Base.x + Obj.B * Base.y +
      Obj.C * Base.z+ Obj.D) / denom;
  return -1;
}

double IntersectSphere(OBJECT Obj, VECTOR Base, VECTOR Dir)
{
  double Det, Aq, Bq, Cq, First, Second;
  VECTOR Base2, Loc;

  Loc.x = Obj.A;  Loc.y = Obj.B;  Loc.z = Obj.C;
  Base2 = Subtract(Base, Loc);
  Aq = Dot(Dir, Dir);
  Bq = 2 * Dot(Dir, Base2);
  Cq = Dot(Base2, Base2) - Obj.D * Obj.D;
  Det = Bq * Bq - 4 * Aq * Cq;
  if (Det >= 0) 
  {
    Det = sqrt(Det);
    First = (-Bq + Det) / (2 * Aq);
    Second = (-Bq - Det) / (2 * Aq);
    if (First < TOL && Second < TOL)
      return -1;
    if (First < Second)
    {
      if (First < TOL)
        return Second;
    }
    else if (Second > TOL)
      return Second;
    return First;
  }
  return -1;
}

// Return the normal N at the point hitPoint for the object.
VECTOR DetermineNormal(OBJECT Obj, VECTOR hitPoint)
{
  VECTOR N;
  switch (Obj.Shape) 
  {
    case PLANE:
      N.x = Obj.A;
      N.y = Obj.B;
      N.z = Obj.C;
      break;
    case SPHERE:
      N.x = hitPoint.x - Obj.A;
      N.y = hitPoint.y - Obj.B;
      N.z = hitPoint.z - Obj.C;
      break;
  }
  N=Normalize(N);
  return N;
}

int InShadow(int hitObject, VECTOR Ray, VECTOR hitPoint)
{
  int i;
  double T;

  for (i=0; i<NUMOBJECTS; i++)
  {  // Check for intersection
    if (i != hitObject)
    {            // Don't check against itself
      switch (Objects[i].Shape)
      {
        case PLANE:
          T = IntersectPlane(Objects[i], hitPoint, Ray);
          break;
        case SPHERE:
          T = IntersectSphere(Objects[i], hitPoint, Ray);
          break;
      }
      if (T > TOL && T < 1)
        return 1;
    }
  }
  return 0;
}

COLOR shade(int Level, int hitObject, VECTOR hitPoint, VECTOR NormalN, VECTOR Dir)
{
  VECTOR L;             
  VECTOR R;
  VECTOR raydir;             
  double NDotL;          
  double DotProd;
  double CosPhi;
  OBJECT Obj;
  COLOR TempI,I;
  
  Obj = Objects[hitObject];

  // Calculate the ambient portion of the point's color
  I.r = Obj.kd.r * Iambient.r;
  I.g = Obj.kd.g * Iambient.g;
  I.b = Obj.kd.b * Iambient.b;

  L.x = Light.position.x - hitPoint.x;   // Calculate the vector from
  L.y = Light.position.y - hitPoint.y;   // the intersection point to
  L.z = Light.position.z - hitPoint.z;   // the light source
  L=Normalize(L);
  
  raydir.x=Dir.x;
  raydir.y=Dir.y;
  raydir.z=Dir.z;
  raydir=Normalize(raydir);


  // Calculate the diffuse portion of the point's color
  if (!InShadow(hitObject, L, hitPoint))
  {
    
    NDotL = Dot(L, NormalN); // Calculate angle between light

    if (NDotL > 0)
	{
      I.r += Obj.kd.r * Light.Color.r * NDotL; 
      I.g += Obj.kd.g * Light.Color.g * NDotL; 
      I.b += Obj.kd.b * Light.Color.b * NDotL;
     }
     
    // light reflection
    R.x = -L.x + 2 * NDotL * NormalN.x;
    R.y = -L.y + 2 * NDotL * NormalN.y;
    R.z = -L.z + 2 * NDotL * NormalN.z;
    R=Normalize(R);

    CosPhi = -Dot(R, raydir);
    if (CosPhi > 0)
	{
      // Add in the specular portion of the point's color
      I.r += Light.Color.r * Obj.ks * pow(CosPhi, Obj.Nzero);
      I.g += Light.Color.g * Obj.ks * pow(CosPhi, Obj.Nzero); 
      I.b += Light.Color.b * Obj.ks * pow(CosPhi, Obj.Nzero); 
    }
  }
 
  if (Level+1 <= MAXLEVEL)
  {
    // recursively trace the ray reflection vector
    DotProd = -Dot(Dir, NormalN);
    R.x = Dir.x + 2 * DotProd * NormalN.x;
    R.y = Dir.y + 2 * DotProd * NormalN.y;
    R.z = Dir.z + 2 * DotProd * NormalN.z;
    R=Normalize(R);

    TempI=TraceRay(Level+1, hitPoint, R);
    I.r += TempI.r * Obj.ks;
    I.g += TempI.g * Obj.ks;
    I.b += TempI.b * Obj.ks;
  }
  return I;
}

COLOR TraceRay(int Level, VECTOR Base, VECTOR Dir)
{
  VECTOR hitPoint, N;

  int i, Closest = -1;
  double NearInt = 100000000000.0, t;

  // Calculate the intersection with the ray and each objec
  for (i=0; i<NUMOBJECTS; i++)
  {
    switch (Objects[i].Shape)
    {
      case PLANE:
        t = IntersectPlane(Objects[i], Base, Dir);
        break;
      case SPHERE:
        t = IntersectSphere(Objects[i], Base, Dir);
        break;
    }
    if (t < NearInt && t > TOL)
    {
      NearInt = t;
      Closest = i;
    }
  }
  
  if (Closest >= 0) 
  { 
    // Calculate the actual intersection point
    hitPoint.x = Base.x + Dir.x * NearInt;
    hitPoint.y = Base.y + Dir.y * NearInt;
    hitPoint.z = Base.z + Dir.z * NearInt;
        
    N=DetermineNormal(Objects[Closest], hitPoint);

    return shade(Level, Closest, hitPoint, N, Dir);
  }
  else
    return background;
}


void myDisplay(void)
{ 
  COLOR I;
  
  VECTOR  Dir;

  VECTOR eye,vup,lookat,w,u,v,scrn;
  double fov,view,scrnu,scrnv,scrnw;
  int ui,vj;

  // set eye position
  eye.x=7.50;
  eye.y=2.5;
  eye.z=10;

  // set lookat position
  lookat.x=0;
  lookat.y=0;
  lookat.z=0;

  // set up vector
  vup.x=0;
  vup.y=1;
  vup.z=0;

  // set field of view
  fov=50;
  // build screen uv coodinate system
  w=Subtract(lookat,eye);
  w=Normalize(w);

  u=Cross(w,vup);
  u=Normalize(u);

  v=Cross(u,w);
  
  // compute the u (and v) coordinate limits for that field of view
  // screen view area will be +/- view in u and v
  view=tan((fov/2)/180*3.14159);

  // Set up viewing coordinate system

  glClear(GL_COLOR_BUFFER_BIT);
  glBegin(GL_POINTS);
  
  // scan all the dots of the "screen"
  for (ui=1;ui<=SCRNWIDTH;ui++)
  {
    for (vj=1;vj<=SCRNHEIGHT;vj++)
	{ //scale sample ui,vj point onto u,v screen point
      scrnu=-view+ui*2*view/SCRNWIDTH;
      scrnv=-view+vj*2*view/SCRNHEIGHT;
      scrnw=1;
      
      //convert screen point to x,y,z
      scrn.x=scrnu*u.x+scrnv*v.x+scrnw*w.x+eye.x;
      scrn.y=scrnu*u.y+scrnv*v.y+scrnw*w.y+eye.y;
      scrn.z=scrnu*u.z+scrnv*v.z+scrnw*w.z+eye.z;

      // build the ray
      Dir=Subtract(scrn,eye);
      
      // fire the ray
      I=TraceRay(1, eye, Dir);
      
      // clamp returned color
      if (I.r > 1.0) I.r = 1.0;
      if (I.g > 1.0) I.g = 1.0;
      if (I.b > 1.0) I.b = 1.0;

      // draw the dot
      glColor3f(I.r,I.g,I.b);
	  glVertex2f(ui,vj);
    }
  }
	 
  glEnd();
  glFlush();
}

int main(int argc, char** argv)
{
/* init the toolkit */
	glutInit(&argc,argv);
/* make a single buffered RGB window */
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB) ;
/* define the initial window size */	
	glutInitWindowSize(SCRNWIDTH,SCRNHEIGHT);
/* define the initial window position */
	glutInitWindowPosition(0,0);
/* title the window */
	glutCreateWindow("Simple Ray Tracing");
/* register the callback function to display */
	glutDisplayFunc(myDisplay);
/* set the background color */
    glClearColor (0.1, 0.2, 0.4, 0.0);
/* set the viewing parameters */
    glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
    gluOrtho2D(0, SCRNWIDTH, 0, SCRNHEIGHT); 
/* go into the perpetual event driven loop */	
	glutMainLoop();
}
