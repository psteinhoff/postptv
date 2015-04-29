int round				(double x);
FILE *fopen_r			(char *filename);
void readseqtrackcrit	(char base_path[256]);

void compose_name_plus_nr				(char *basename,char *str,int nr,char *filename);
void compose_name_plus_nr_str			(char *basename,char *str,int nr,char *filename);
void showErrorMessageOnFileOpening		(char *filename);
int countConnectedParticles				(int size1, vector *line1, int size2, vector *line2);
int countValidatedParticleConnections	(vector *line, int upperLimit);

/*functions to print to cmd*/
void printSimpleText	(int x, int y, char const *text);
void printSimpleLine	(int endOfline,int line,int textOffset, char const *newText);
void printHeader		();
void progressBar		(int line, int width, int progressInPercent);
void clearLine			(int line, int width);
