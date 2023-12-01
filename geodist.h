/* Biblioteca com funções para cálculo de distâncias
 * */
#ifndef _GEODIST_
#define _GEODIST_


#define F_ERROR         (0.0001)                    //Erro permitido (float)
#define EQ_FLOAT(a,b)   (fabs((a)-(b)) < F_ERROR)   //Se dois floats são iguais

#define BETA            2       //Potência utilizada na interpolação (valor maior deixa mais suave)
#define MAX_DIST        3       //Distância máxima de busca por quadricula (em latlon)

#define EARTH_RADIUS    6371.0  //Radius of the Earth in kilometers


/* Conversão de graus (latitude e longitude) para radianos*/
double to_radians(double degree);

/* Retorna a distância euclidiana entre A e B, sendo
 * A = (ax,ay) e B = (bx,by)
**/
double distance(double ax, double ay, double bx, double by);

/* Retorna a distância euclidiana ao quadrado entre A e B, sendo
 * A = (ax,ay) e B = (bx,by)
**/
double sqr_distance(double ax, double ay, double bx, double by);

/* Retorna a distância na esfera entre os pares de (lat,lon) A e B
 * Haversine Distance (https://en.wikipedia.org/wiki/Haversine_formula)
 * A = (ax,ay) e B = (bx,by)
**/
double haversine_distance(double loni, double lati, double lonf, double latf);

/* Função auxiliar para verificar se uma coordenada está dentro da área 
 * Retorna 1 se a coordenada (x,y) está dentro da área delimitada por (xi,xf,yi,yf).
 * Retorna 0 caso contrário.
**/
int inside_area(double x, double y, double xi, double xf, double yi, double yf);

/* Função auxiliar para verificar se um valor está dentro do intervalo
 * Retorna 1 se o valor 'val' está contido no intervalo (inclusivo).
 * Retorna 0 caso contrário.
**/
int inside_axis(double val, double initial, double final);



#endif
