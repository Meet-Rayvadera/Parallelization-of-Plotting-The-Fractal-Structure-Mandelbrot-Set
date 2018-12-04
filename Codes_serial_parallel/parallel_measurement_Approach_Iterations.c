/*
	Utsav Patel: 201501443
	Meet Rayvadera: 201501444
	"Mandelbrot Set"
*/

/*
	Parallel code for dividing iterations across different cores (Pipelining approach through Queue)
*/

#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<time.h>
#include<string.h>
#include<stdlib.h>

struct QNode
{
    double x,y;
	int finish;
    struct QNode *next;
};
 
// The queue, front stores the front node of LL and rear stores the
// last node of LL
struct Queue
{
    struct QNode *front, *rear;
} *qu[20];
 
// A utility function to create a new linked list node.
struct QNode* newNode(double a,double b)
{
    struct QNode *temp = (struct QNode*)malloc(sizeof(struct QNode));
    temp->x = a;
	temp->y = b;
	temp->finish = 0;
    temp->next = NULL;
    return temp; 
}
 
// A utility function to create an empty queue
struct Queue *createQueue()
{
    struct Queue *q = (struct Queue*)malloc(sizeof(struct Queue));
    q->front = q->rear = NULL;
    return q;
}	

void enQueue(struct Queue *q,struct QNode *temp)
{
    // Create a new LL node
 
    // If queue is empty, then new node is front and rear both
    if (q->rear == NULL)
    {
       q->front = q->rear = temp;
       return;
    }
 
    // Add the new node at the end of queue and change rear
    q->rear->next = temp;
    q->rear = temp;
}
 
// Function to remove a key from given queue q
struct QNode *deQueue(struct Queue *q)
{
    // If queue is empty, return NULL.
    if (q->front == NULL)
       return NULL;
 
    // Store previous front and move front one node ahead
    struct QNode *temp = q->front;
    q->front = q->front->next;
 
    // If front becomes NULL, then change rear also as NULL
    if (q->front == NULL)
       q->rear = NULL;
    return temp;
}

int isEmpty(struct Queue *q){
	if(q->front==NULL)
		return 1;
	return 0;
}

typedef struct {
  unsigned char red,green,blue;
} PPMPixel;

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

	radius = 2;
	
	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
	imageWidth = n;
	imageHeight = n;
	minX = -2;
	maxX = 1;
	minY = -1.5;
	maxY = 1.5;
	maxIteration = 120;
	
	char *problem_name = "mandelbrot_set";
	char *approach_name = "Parallelism of maxIteration";
	
	FILE* outputFile;
	char outputFileName[50];		
	sprintf(outputFileName,"output/%s_%s_%s_%s_output.txt",problem_name,approach_name,argv[1],argv[2]);
	
	long int Iterations[imageHeight][imageWidth];
	int i,j;
	
	PPMImage* output = (PPMImage *)malloc(sizeof(PPMImage));
	output->x = imageHeight;
	output->y = imageWidth;
	output->data = (PPMPixel *) malloc(imageHeight*imageWidth*sizeof(PPMPixel));	

	clock_gettime(CLK, &start_alg);	/* Start the algo timer */

int numth = omp_get_num_threads(); // 7
int seq = maxIteration/numth; // 120/7 = 17

if(maxIteration%numth != 0)  // 120%7 = 1
	seq++;  // 7++ = 8

for(i=0;i<20;i++)
	qu[i] = createQueue();

#pragma omp parallel private(i,j) shared(qu,seq,numth)
{
		int yi = omp_get_thread_num();
		for(i=0; i<imageHeight; i++){
			for(j=0; j<imageWidth; j++){
			
			
			if(yi!=0 && isEmpty(qu[yi])){    // run infinite loop until we get some pre calculate values
				j--;
				continue;
			}
	        printf("%d %d %d\n",yi,i,j);
			double x = 0,y = 0;
			if(yi!=0){
				struct QNode *aa = deQueue(qu[yi]);          //  if there exist some pre calculated value, we get these value from queue
				if(aa){
					if(aa->finish){
						struct QNode *bb = newNode(x,y);            // if x and y are already out of bound then this finish variable is true and this thred don't run the loop
						bb->finish = 1;
						enQueue(qu[yi+1],bb);
						continue;
					}
					x = aa->x;
					y = aa->y;
				}
			}
			int index = i*imageWidth + j;
			PPMPixel *tem = output->data + index;
			
			double xtemp;
			
			int start = yi*seq;
			int end = start+seq;
			if(end>maxIteration)
				end = maxIteration;
			
			int k = 0;
			int cc = 0;
			for(k=start;k<end;k++){                                   // take precomputed x and y values and run loop for this thread and find values
				xtemp = x*x - y*y + convertedX(j);
				y = 2*x*y + convertedY(i);
				x = xtemp;
		
				//printf("%f + i*%f,    %ld\n",x,y,t);
		
				if(x*x + y*y > radius*radius){ // i.e. If  point is not in circle of radius 2(or any other value) centred at 0
					Iterations[i][j] = k;
					struct QNode *bb = newNode(0,0);
					bb->finish = 1;
					enQueue(qu[yi+1],bb);
					cc = 1;
					break; 
				}
			}
			
			if(!cc){
				struct QNode *bb = newNode(x,y);
				enQueue(qu[yi+1],bb);
			}
			
			else{
				long int tem2 = Iterations[i][j];
				tem2 = 255.0 - ((tem2*255.0)/(double) maxIteration);  //Normalize number of iterations from range 1 to maxiterations to range 0 to 255
				tem->red = tem2;                            //All 3 are random choice by try and errror
				tem->green = tem2*8;
				tem->blue = tem2*4;
			}
		}
	}
}

	clock_gettime(CLK, &end_alg);	/* End the algo timer */

	writePPM("poutput.ppm", output);

	clock_gettime(CLK, &end_e2e);

	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);
	
	double e = e2e.tv_sec + (1e-9)*e2e.tv_nsec;
	double a = alg.tv_sec + (1e-9)*alg.tv_nsec;

	// printf("%d,  %ld ,  %d ,  %ld \n", start_alg.tv_sec, start_alg.tv_nsec, end_alg.tv_sec, end_alg.tv_nsec);
	// printf("%d,  %ld ,  %d ,  %ld \n %.10f %.10f \n", e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec, e, a);
	
	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, n, p, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);
	return 0;
}