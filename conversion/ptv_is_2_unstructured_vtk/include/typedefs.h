#define sqr(x) x*x
#define maxcand 60//Beat changed it on 090325 // P.Steinhoff 121211 handle with care!! Can cause problems with allocation
#define maxtrack 64 //Beat changed it on 090325

#ifndef Ms
#define Ms 5000 //Beat changed it on 090325 // P.Steinhoff 120302 changed from 20000 to 40000 
#endif

#ifndef  M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define POSI 256 //Beat changed it on 090325

typedef	double	Dmatrix[3][3];	/* 3 x 3 rotation matrix */

typedef struct
{
  int pnr;
  double x, y, z;
    }
coord_3d;

typedef struct
{
  int pnr;
  double x, y;
}
coord_2d;

typedef struct
{
  double  x0, y0, z0;
  double  omega, phi, kappa;
  Dmatrix dm;
}
Exterior;

typedef struct
{
  double  dacc, dangle, dvxmax, dvxmin;
  double dvymax, dvymin, dvzmax, dvzmin;
  int  dsumg, dn, dnx, dny, add;
  int showAdds;
  int nframe, deltaT;
  float link_tol, jump_tol, Vmax;
  int flag[8];
  int glue,track,iter,remKints,writeBin,writeTraj,out,sameCamConf;
}
trackparameters;

typedef struct
{
  double xh, yh;
  double cc;
}
Interior;

typedef struct
{
  double vec_x,vec_y,vec_z;
}
Glass;

typedef struct
{
  double k1,k2,k3,p1,p2,scx,she;
}
ap_52;

typedef struct
{
  double n1,n2,n3,d;
  int    lut;
}
mm_3p;

typedef struct
{
  int  	  nlay;
  double  n1;
  double  n2[3];
  double  d[3];
  double  n3;
  int     lut;
}
mm_np;

typedef struct
{
	int j;
	double nx_ratio;
	double ny_ratio;
	double n_ratio;
	double sumg_ratio;
	double nxy_ratio;
	double average_gv_ratio;
}
spatio_temp_particle_optix;

typedef struct
{
	double dl;
	double acc;
	double Acc[2];
	double angle;
	double max_angle;
	double sumg[2];
	double n[2];
	double nx[2];
	double ny[2];
	double coefficients[6];
	int prev;
}
spring_parameters;

typedef struct
{
  int     pnr;
  double  x, y;
  int     n, nx, ny, sumg;
  int     tnr;
}
target;

typedef struct
{
  int 	pos, status;
  short	xmin, xmax, ymin, ymax;
  int   n, sumg;
  double  x, y;
  int   unr, touch[4], n_touch;	/* unified with target unr, touching ... */
}
peak;

typedef struct
{
  short	       	x,y;
  unsigned short	g; //changed to unsigned short due to work with 16 bit values
  short	       	tnr;
}
targpix;

typedef struct
{
  int  	left, right;
  double  tol, corr;
}
conjugate;

typedef struct
{
  int  	pnr;
  double  tol, corr;
}
candidate;

typedef struct
{
  int     p[4];
  double  corr;
}
n_tupel;

typedef struct
{
  int nr;
  int p[4];
}
corres;

typedef struct
{
  int    	p1;	       	/* point number of master point */
  int    	n;	       	/* # of candidates */
  int    	p2[maxcand];	/* point numbers of candidates */
  double	corr[maxcand];	/* feature based correlation coefficient */
  double	dist[maxcand];	/* distance perpendicular to epipolar line */
}
correspond;	       	/* correspondence candidates */

typedef struct
{
  coord_3d	origin;
  int          	nr,nz;
  double       	rw;
  double       	*data;
}
mm_LUT;


typedef struct /* struct for what was found to corres */
{
 int ftnr, freq, whichcam[4];
}
foundpix;

typedef struct
{
  int multi, h[maxtrack], freq[maxtrack];
  double quali[maxtrack];
}
currlist;

typedef struct
{
  int p,n;
  double x1, y1, z1;
  double Ux, Uy, Uz, Um;
  int type;
}
vector;

typedef struct
{
  int z1, c[maxcand], n[maxcand];
  double quali[maxcand];
}
prevlist;

typedef struct Pstruct
{
  float x[3]; /*coordinates*/
  int prev, next; /*pointer to prev or next link*/
  int prio; /*Prority of link is used for different levels*/
  float decis[POSI]; /*Bin for decision critera of possible links to next dataset*/
  float finaldecis; /*final decision critera by which the link was established*/
  int linkdecis[POSI]; /* pointer of possible links to next data set*/
  int inlist; /* Counter of number of possible links to next data set*/
  float EBC_decis[POSI];
} P;

typedef struct trackingAttributes
{
	int sumgMin[4];
	int sumgMax[4];
}A;

typedef struct img16uStruct
{ //Similar to OpenCV
	unsigned short *pixelData;
	int PixValue[3000][3000];
	int img_width;
	int img_height;
	int depth;		/* Pixel depth in bits*/
	int nChannels;
	int imageSize;  /* Image data size in bytes (==image->height*image->widthStep
                               in case of interleaved data)*/
	int widthStep;	/* Size of aligned image row in bytes.    */
	//img16uStruct Img16u;
	//Img16u = (unsigned short*) malloc((unsigned short)(fileImgHeight*fileImgWidth);
} Img16U;

#define MAX_FILENAME_LEN 1024
#define FILENAME_IN		"res/rt_is"
#define FILENAME_OUT	"res/ptv_is"
#define PARAMETERFILE	"res/track.par"
#define COORDPARAFILE	"res/coord.par"
#define STATUSFILE		"res/track.out"

#define velo_par		"parameters/velo_trans.par"
#define sequence_par	"parameters/sequence.par"
#define track_par		"parameters/track.par"
#define criteria_par	"parameters/criteria.par"
#define point_offset	"parameters/3d_point_offset.par"
