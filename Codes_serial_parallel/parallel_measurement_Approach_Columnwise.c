/*
	Utsav Patel: 201501443
	Meet Rayvadera: 201501444
	"Mandelbrot Set"
*/

/*
	Parallel Code for column-wise parallelization
*/

#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<time.h>
#include<string.h>
#include<stdlib.h>

//Structure for each pixel
typedef struct {
  unsigned char red,green,blue;
} PPMPixel;

//Structure for PPM image
typedef struct {
  int x, y;
  PPMPixel *data;
} PPMImage;

#define RGB_COMPONENT_COLOR 255

//  Using the MONOTONIC clock 
#define CLK CLOCK_MONOTONIC

/* Function to compute the difference between two points in time */
struct timespec diff(struct timespec start, struct timespec end);

struct timespec diff(struct timespec start, struct timespec end)
{
	struct timespec temp;
	if((end.tv_nsec - start.tv_nsec) < 0){
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}
	return temp;
}

//For reading PPM image
static PPMImage *readPPM(const char *filename)
{
  char buff[16];
  PPMImage *img;
  FILE *fp;
  int c, rgb_comp_color;
  //open PPM file for reading
  fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    exit(1);
  }

  //read image format
  if (!fgets(buff, sizeof(buff), fp)) {
    perror(filename);
    exit(1);
  }

  //check the image format
  if (buff[0] != 'P' || buff[1] != '6') {
    fprintf(stderr, "Invalid image format (must be 'P6')\n");
    exit(1);
  }

  //alloc memory form image
  img = (PPMImage *)malloc(sizeof(PPMImage));
  if (!img) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  //check for comments
  c = getc(fp);
  while (c == '#') {
    while (getc(fp) != '\n') ;
    c = getc(fp);
  }

  ungetc(c, fp);
  //read image size information
  if (fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
    fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
    exit(1);
  }

  //read rgb component
  if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
    fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
    exit(1);
  }

  //check rgb component depth
  if (rgb_comp_color!= RGB_COMPONENT_COLOR) {
    fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
    exit(1);
  }

  while (fgetc(fp) != '\n') ;
  //memory allocation for pixel data
  img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));

  if (!img) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  //read pixel data from file
  if (fread(img->data, 3 * img->x, img->y, fp) != img->y) {
    fprintf(stderr, "Error loading image '%s'\n", filename);
    exit(1);
  }

  fclose(fp);
  return img;
}

//For writing into PPM image
void writePPM(const char *filename, PPMImage *img)
{
  FILE *fp;
  //open file for output
  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    exit(1);
  }

  //write the header file
  //image format
  fprintf(fp, "P6\n");

  //comments


  //image size
  fprintf(fp, "%d %d\n",img->x,img->y);

  // rgb component depth
  fprintf(fp, "%d\n",255);

  // pixel data
  fwrite(img->data, 3 * img->x, img->y, fp);
  fclose(fp);
}




double radius;  // Bounding value for absolute value of terms of sequence for each point     
int imageHeight,imageWidth; // Number of pixels of image
long int maxIteration;  // Number of terms in sequence for each point
double minX,maxX,minY,maxY; // Range of points in complex plane


double convertedX(int x){			// Normalization of point ranging between 0 to imageWidth-1 to range between minX and maxX  
	if(imageWidth == 1)		return (minX + maxX)/2;
	
	double dis = maxX-minX;
	double d = minX + (x*dis)/(imageWidth-1);
	return d;
}

double convertedY(int y){		// Normalization of point ranging between 0 to imageHeight-1 to range between minY and maxY  
	if(imageHeight == 1)		return (minY + maxY)/2;
	
	double dis = maxY-minY;
	double d = maxY - (y*dis)/(imageHeight-1);
	return d;
}

long int findIterations(double a,double b){ // This function return the iterations till our sequence remains bounded.    
	long int t = 0;
	double x = 0, y = 0, xtemp;
	
	while(++t < maxIteration){
		xtemp = x*x - y*y + a; // xnew + i*ynew = (x*x - y*y + a) + i*(2*x*y + b)
		y = 2*x*y + b;
		x = xtemp;
		
		if(x*x + y*y > radius*radius) // i.e. If  point is not in circle of radius 2(or any other value) centred at 0
			break;                    
	}
	return t;  
}

int main(int argc, char* argv[]){

	struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
	/* Should start before anything else */
	clock_gettime(CLK, &start_e2e);	

	int n = atoi(argv[1]);	
	int p = atoi(argv[2]);
	imageWidth = n;
	imageHeight = n;
	minX = -2;
	maxX = 1;
	minY = -1.5;
	maxY = 1.5;
	maxIteration = 120;
	radius = 2;

	char *problem_name = "mandelbrot_set";
	char *approach_name = "ColumnWise";

	FILE* outputFile;
	char outputFileName[50];		
	sprintf(outputFileName,"output/%s_%s_%s_%s_output.txt",problem_name,approach_name,argv[1],argv[2]);	
	
	long int Iterations[imageHeight][imageWidth];
	int i,j;
	
	PPMImage* output = (PPMImage *)malloc(sizeof(PPMImage));
	output->x = imageHeight;
	output->y = imageWidth;
	output->data = (PPMPixel *) malloc(imageHeight*imageWidth*sizeof(PPMPixel));	

	omp_set_num_threads(p); // Sets number of threads
	clock_gettime(CLK, &start_alg);	/* Start the algo timer */

#pragma omp parallel private(i,j) shared(Iterations,output)
{
#pragma omp for // Divide columns across different cores
	for(j=0; j<imageWidth; j++){
		for(i=0; i<imageHeight; i++){
			
			int index = i*imageWidth + j;
			PPMPixel *tem = output->data + index;
			
			Iterations[i][j] = findIterations(convertedX(j),convertedY(i));  //Writing into shared memory 1
			long int tem2 = Iterations[i][j];
			
			tem2 = 255.0 - ((tem2*255.0)/(double) maxIteration);  //Normalize number of iterations from range 1 to maxiterations to range 0 to 255
			tem->red = tem2;                            //All 3 are random choice by try and errror     //Writing into shared memory 2
			tem->green = tem2*8;						//Writing into shared memory 3
			tem->blue = tem2*4;							//Writing into shared memory 4
			
		}
	}
}

	clock_gettime(CLK, &end_alg);	/* End the algo timer */

	writePPM("poutput.ppm", output);

	clock_gettime(CLK, &end_e2e);

	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);
	
	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, n, p, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;
}
