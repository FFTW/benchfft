/* MODIFIED: 8/17/98 by Steven G. Johnson (stevenj@alum.mit.edu) for
   inclusion in benchfft (http://www.fftw.org/benchfft).  I
   changed data types to match those in the benchmark, and changed
   some subroutine names to prevent conflicts.  Modified 9/5/98 to
   use ANSI C function declarations instead of K&R. */

/* MODIFIED 8/31/02 by Steven G. Johnson to use header and data types
   from nbenchfft.  The original site and code seem to have disappeared
   from the web.  See also NOTE below on accuracy. */

#include "bench-user.h"
typedef bench_real real;

#define MY_PI 3.1415926535897932384626434

/************************************************************************/
/*	fft generalisee: Yves MONNIER 1995	                        */
/*      http://www.igd.u-bordeaux.fr/~monnier/fft.html                  */ 
/*                                                                      */
/*	compilation sous unix:	cc fft.c -o fft -lm               	*/
/*	Pour vous en servir, rien de plus simple. D'abord faire un	*/
/*	appel a la fonction ouvre_fft() en lui donnant le nombre	*/
/*	d'elements contenu dans le signal. Ensuite appeler autant	*/
/*	de fois que desire la fonction fft() avec en parametre le	*/
/*	tableau contenant le signal (suite de valeurs complexes		*/
/*	definies en float: reel,imaginaire,reel,imaginaire,reel,...)	*/
/*	Le tableau contient en sortie le spectre du signal.		*/
/*	La derniere transformee faite, appelez la ferme_fft()		*/
/*	UN MAIN EST FOURNI PLUS BAS QUI DONNE LA DEMARCHE A SUIVRE	*/
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIRECTE		 1.0
#define INVERSE		-1.0
#define OUVERT		 1
#define FERME		 0

#define W1_3r 		-0.5
#define W1_3i 		 0.8660254037844386467637231707529361834715
#define W1_5r 		 0.3090169943749474241022934171828190588603
#define W1_5i 		 0.951056516295153572116439333378
#define W2_5r 		-0.80901699437494742410229341718281
#define W2_5i 		 0.5877852522924731291687059546391

static int facteurs_premier[]={97,89,83,79,73,71,67,61,59,53,47,43,41,37,31,29,23,19,17,13,11,7,5,4,3,2,0}; /* ne pas mettre 1 */

static char	 fft_flag=FERME;
static real	 fft_sens;
static int 	 fft_taille;
static real 	*fft_tampon;
static int 	*fft_facteurs;
static int 	*fft_poids;
static int	*fft_arrangement;
static real	*fft_wn;

/******************************************************************************/

int fabrique_poids_et_facteurs(void)
{
int n;
int i,j;
int nb_facteurs;

n=fft_taille;
nb_facteurs=0;

for(i=0;facteurs_premier[i];i++)
	{
	while(n%facteurs_premier[i]==0)
		{
		n/=facteurs_premier[i];
		nb_facteurs++;
		}
	}
if(n!=1)
	nb_facteurs ++;

if((fft_facteurs=(int*)malloc((nb_facteurs+1)*sizeof(int)))==NULL)
	return (-1);
if((fft_poids=(int*)malloc((nb_facteurs+1)*sizeof(int)))==NULL)
	{
	free(fft_facteurs);
	return(-1);
	}
j=0;
fft_poids[0]=fft_taille;
if(n!=1)
	{
	fft_facteurs[0]=n;
	fft_poids[1]=fft_taille/n;
	j=1;
	}
n=fft_taille;
for(i=0;facteurs_premier[i];i++)
	{
	while(n%facteurs_premier[i]==0)
		{
		n/=(fft_facteurs[j++]=facteurs_premier[i]);
		fft_poids[j]=fft_poids[j-1]/facteurs_premier[i];
		}
	}
fft_facteurs[j]=0;
return(0);
}

/******************************************************************************/

int fabrique_le_tableau_arrangement(void)
{
int n,i,j;

if((fft_arrangement=(int*)malloc(fft_taille*sizeof(int)))==NULL)
	return(-1);
n=0;
for(i=0;i<fft_taille;i++)
	{
	fft_arrangement[i]=2*n;
	n+=fft_poids[1];
	j=1;
	while(n>=fft_poids[j-1])
		{
		n+=fft_poids[j+1]-fft_poids[j-1];
		j++;
		}
	}
return(0);
}

/******************************************************************************/

void arrange(real *donnees)
{
int i,j;
real f;

for(i=0,j=0;i<fft_taille;i++,j+=2)
	{
	fft_tampon[j]=donnees[fft_arrangement[i]];
	fft_tampon[j+1]=donnees[fft_arrangement[i]+1];	
	}
for(i=0;i<fft_taille*2;i++)
	donnees[i]=fft_tampon[i];
}

/******************************************************************************/

int fabrique_le_tableau_wn(void)
{
int i;
real n,m;

n=fft_sens*2.0*MY_PI/(real)fft_taille;
m=0.0;
if((fft_wn=(real*)malloc(fft_taille*2*sizeof(real)))==NULL)
	return(-1);
for(i=0;i<fft_taille*2;i+=2)
	{
	     /* NOTE by Steven G. Johnson, 8/31/02: computing the
		phase angle by repeated additions here completely
		kills the accuracy.  (e.g. in single precision the FFT
		accuracy for 2^16 is only ~2e-3!)  Using i*(0.5*n) is
		just as fast, and is way more accurate. (e.g. the
		abovementioned 2^16 FFT accuracy becomes ~3e-7.) 
	        We preserve the author's version as a cautionary tale. */
	fft_wn[i]=cos(m);
	fft_wn[i+1]=sin(m);
	m+=n;
	}
return(0);
}

/******************************************************************************/

int ouvre_fft(int taille,real sens)
{
if(fft_flag==OUVERT)
	return(-1);
fft_taille=taille;
if((fft_tampon=(real *)malloc(taille*2*sizeof(real)))==NULL)
	return(-1);
if(fabrique_poids_et_facteurs()==-1)
	{
	free(fft_tampon);
	return(-1);
	}
if(fabrique_le_tableau_arrangement()==-1)
	{
	free(fft_tampon);
	free(fft_facteurs);
	free(fft_poids);
	return(-1);
	}
fft_sens=sens;
if(fabrique_le_tableau_wn()==-1)
	{
	free(fft_tampon);
	free(fft_arrangement);
	free(fft_facteurs);
	free(fft_poids);
	return(-1);
	}
fft_flag=OUVERT;
return(0);
}

/******************************************************************************/

int ferme_fft(void)
{
if(fft_flag==FERME)
	return(-1);
free(fft_tampon);
free(fft_arrangement);
free(fft_facteurs);
free(fft_poids);
free(fft_wn);
fft_flag=FERME;
return(0);
}

/******************************************************************************/

int fft(real *donnees)
{
int i,entree,indice;
int paquet,module,nb_modules;
int er0,ei0,er1,ei1,er2,ei2,er3,ei3,er4,ei4;
int ern,ein;
int wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4;
int w,wr0,wrn,win;
real dr0,di0,dr1,di1,dr2,di2,dr3,di3,dr4,di4;
real frn,fin;
real Ar,Ai,Br,Bi,Cr,Ci,Dr,Di;
real f;

if(fft_flag==FERME)
	return(-1);

arrange(donnees);

for(i=0;fft_facteurs[i];i++)
	{
	for(paquet=0;paquet<2*fft_taille;paquet+=2*fft_taille/fft_poids[i+1])
		{
		nb_modules=2*fft_taille/fft_poids[i];
		for(module=0;module<nb_modules;module+=2)
			{
			switch(fft_facteurs[i])
				{
				case 5 :	er0=module+paquet;
						ei0=er0+1;
						er1=er0+nb_modules;
						ei1=er1+1;
						er2=er1+nb_modules;
						ei2=er2+1;
						er3=er2+nb_modules;
						ei3=er3+1;
						er4=er3+nb_modules;
						ei4=er4+1;
						wr1=module*fft_poids[i+1];
						wi1=wr1+1;
						wr2=wr1<<1;
						wi2=wr2+1;
						wr3=wr1*3;
						wi3=wr3+1;
						wr4=wr1<<2;
						wi4=wr4+1;
						dr0=donnees[er0];
						di0=donnees[ei0];
						dr1=donnees[er1]*fft_wn[wr1]-donnees[ei1]*fft_wn[wi1];
						di1=donnees[er1]*fft_wn[wi1]+donnees[ei1]*fft_wn[wr1];
						dr2=donnees[er2]*fft_wn[wr2]-donnees[ei2]*fft_wn[wi2];
						di2=donnees[er2]*fft_wn[wi2]+donnees[ei2]*fft_wn[wr2];
						dr3=donnees[er3]*fft_wn[wr3]-donnees[ei3]*fft_wn[wi3];
						di3=donnees[er3]*fft_wn[wi3]+donnees[ei3]*fft_wn[wr3];
						dr4=donnees[er4]*fft_wn[wr4]-donnees[ei4]*fft_wn[wi4];
						di4=donnees[er4]*fft_wn[wi4]+donnees[ei4]*fft_wn[wr4];
						Ar=dr0+W1_5r*(dr1+dr4)+W2_5r*(dr2+dr3);
						Ai=di0+W1_5r*(di1+di4)+W2_5r*(di2+di3);
						Br=W1_5i*(di1-di4)+W2_5i*(di2-di3);
						Bi=W1_5i*(dr1-dr4)+W2_5i*(dr2-dr3);
						Cr=dr0+W2_5r*(dr1+dr4)+W1_5r*(dr2+dr3);
						Ci=di0+W2_5r*(di1+di4)+W1_5r*(di2+di3);
						Dr=W2_5i*(di1-di4)-W1_5i*(di2-di3);
						Di=W2_5i*(dr1-dr4)-W1_5i*(dr2-dr3);
						donnees[er0]=dr0+dr1+dr2+dr3+dr4;
						donnees[ei0]=di0+di1+di2+di3+di4;
						if(fft_sens==DIRECTE)
							{
							donnees[er1]=Ar-Br;
							donnees[ei1]=Ai+Bi;
							donnees[er2]=Cr-Dr;
							donnees[ei2]=Ci+Di;
							donnees[er3]=Cr+Dr;
							donnees[ei3]=Ci-Di;
							donnees[er4]=Ar+Br;
							donnees[ei4]=Ai-Bi;
							}
						else
							{
							donnees[er4]=Ar-Br;
							donnees[ei4]=Ai+Bi;
							donnees[er3]=Cr-Dr;
							donnees[ei3]=Ci+Di;
							donnees[er2]=Cr+Dr;
							donnees[ei2]=Ci-Di;
							donnees[er1]=Ar+Br;
							donnees[ei1]=Ai-Bi;
							}
					break;
				case 4 :	er0=module+paquet;
						ei0=er0+1;
						er1=er0+nb_modules;
						ei1=er1+1;
						er2=er1+nb_modules;
						ei2=er2+1;
						er3=er2+nb_modules;
						ei3=er3+1;
						wr1=module*fft_poids[i+1];
						wi1=wr1+1;
						wr2=wr1<<1;
						wi2=wr2+1;
						wr3=wr1*3;
						wi3=wr3+1;
						dr0=donnees[er0];
						di0=donnees[ei0];
						dr1=donnees[er1]*fft_wn[wr1]-donnees[ei1]*fft_wn[wi1];
						di1=donnees[er1]*fft_wn[wi1]+donnees[ei1]*fft_wn[wr1];
						dr2=donnees[er2]*fft_wn[wr2]-donnees[ei2]*fft_wn[wi2];
						di2=donnees[er2]*fft_wn[wi2]+donnees[ei2]*fft_wn[wr2];
						dr3=donnees[er3]*fft_wn[wr3]-donnees[ei3]*fft_wn[wi3];
						di3=donnees[er3]*fft_wn[wi3]+donnees[ei3]*fft_wn[wr3];
						Ar=dr0+dr2;
						Ai=di0+di2;
						Br=dr0-dr2;
						Bi=di0-di2;
						Cr=dr1+dr3;
						Ci=di1+di3;
						Dr=di1-di3;
						Di=dr1-dr3;
						donnees[er0]=Ar+Cr;
						donnees[ei0]=Ai+Ci;
						donnees[er2]=Ar-Cr;
						donnees[ei2]=Ai-Ci;
						if(fft_sens==DIRECTE)
							{
							donnees[er1]=Br-Dr;
							donnees[ei1]=Bi+Di;
							donnees[er3]=Br+Dr;
							donnees[ei3]=Bi-Di;
							}
						else
							{
							donnees[er3]=Br-Dr;
							donnees[ei3]=Bi+Di;
							donnees[er1]=Br+Dr;
							donnees[ei1]=Bi-Di;
							}
					break;
				case 3 :	er0=module+paquet;
						ei0=er0+1;
						er1=er0+nb_modules;
						ei1=er1+1;
						er2=er1+nb_modules;
						ei2=er2+1;
						wr1=module*fft_poids[i+1];
						wi1=wr1+1;
						wr2=wr1<<1;
						wi2=wr2+1;
						dr0=donnees[er0];
						di0=donnees[ei0];
						dr1=donnees[er1]*fft_wn[wr1]-donnees[ei1]*fft_wn[wi1];
						di1=donnees[er1]*fft_wn[wi1]+donnees[ei1]*fft_wn[wr1];
						dr2=donnees[er2]*fft_wn[wr2]-donnees[ei2]*fft_wn[wi2];
						di2=donnees[er2]*fft_wn[wi2]+donnees[ei2]*fft_wn[wr2];
						Ar=dr0+W1_3r*(dr1+dr2);
						Br=W1_3i*(di2-di1);
						Ai=di0+W1_3r*(di1+di2);
						Bi=W1_3i*(dr1-dr2);
						donnees[er0]=dr0+dr1+dr2;
						donnees[ei0]=di0+di1+di2;
						if(fft_sens==DIRECTE)
							{
							donnees[er1]=Ar+Br;
							donnees[ei1]=Ai+Bi;
							donnees[er2]=Ar-Br;
							donnees[ei2]=Ai-Bi;
							}
						else
							{
							donnees[er2]=Ar+Br;
							donnees[ei2]=Ai+Bi;
							donnees[er1]=Ar-Br;
							donnees[ei1]=Ai-Bi;
							}
					break;
				case 2 :	er0=module+paquet;
						ei0=er0+1;
						er1=er0+nb_modules;
						ei1=er1+1;
						wr1=module*fft_poids[i+1];
						wi1=wr1+1;
						dr0=donnees[er0];
						di0=donnees[ei0];
						dr1=donnees[er1]*fft_wn[wr1]-donnees[ei1]*fft_wn[wi1];
						di1=donnees[er1]*fft_wn[wi1]+donnees[ei1]*fft_wn[wr1];
						donnees[er0]=dr0+dr1;
						donnees[ei0]=di0+di1;
						donnees[er1]=dr0-dr1;
						donnees[ei1]=di0-di1;
					break;
				default : 	ern=module+paquet;
						wr0=module*fft_poids[i+1];
						for(entree=0;entree<fft_facteurs[i]*2;entree+=2)
							{
							ein=ern+1;
							wrn=(wr0*entree/2)%(2*fft_taille);
							win=wrn+1;
							fft_tampon[entree]=donnees[ern]*fft_wn[wrn]-donnees[ein]*fft_wn[win];
							fft_tampon[entree+1]=donnees[ern]*fft_wn[win]+donnees[ein]*fft_wn[wrn];
							ern+=nb_modules;
							}
						w=2*fft_taille/fft_facteurs[i];
						ern=module+paquet;
						for(entree=0;entree<fft_facteurs[i]*2;entree+=2)
							{
							frn=0.0;
							fin=0.0;
							for(indice=0;indice<fft_facteurs[i]*2;indice+=2)
								{
								wrn=(w*entree*indice/4)%(2*fft_taille);
								frn+=fft_tampon[indice]*fft_wn[wrn]-fft_tampon[indice+1]*fft_wn[wrn+1];
								fin+=fft_tampon[indice]*fft_wn[wrn+1]+fft_tampon[indice+1]*fft_wn[wrn];
								}
							donnees[ern]=frn;
							donnees[ern+1]=fin;
							ern+=nb_modules;
							}
					break;
				}
			}
		}
	}
if(fft_sens==INVERSE)
	for(f=(real)fft_taille,i=0;i<fft_taille*2;i++)
		donnees[i]/=f;
return(0);
}

#if 0
/*********Exemple d'utilisation de la fft**********/

/******************************************************************************/

real module(reel,imag,n)
real reel,imag;
int n;
{
return(sqrt((reel*reel)+(imag*imag))/(real)n);
}

int main(argc,argv)
int argc;
char *argv[];
{
int i,j,taille;
real *donnees,frequence,f;

if(argc!=3 && argc!=4)
	{
	printf("usage: %s <nb_echantillons> <frequence> [p]\n",argv[0]);
	return(0);
	}

taille=atoi(argv[1]);
frequence=atof(argv[2]);

if((donnees=(real*)malloc(taille*2*sizeof(real)))==NULL)
	return(-1);

for(i=0;i<taille*2;i+=2)
	{
	donnees[i]=cos((real)i*2*MY_PI/frequence);
	donnees[i+1]=0.0;
	}

if(argc==4)
	{
	for(f=1;f>0;f-=0.05)
		{
		for(i=0;i<taille*2;i+=2)
			if(donnees[i]>f)
				printf("*");
			else
				printf(" ");
		printf("\n");
		}
	for(f=0;f>-1;f-=0.05)
		{
		for(i=0;i<taille*2;i+=2)
			if(donnees[i]<f)
				printf("*");
			else
				printf(" ");
		printf("\n");
		}
	}
printf("SIGNAL D'ENTREE\n");

if(ouvre_fft(taille,DIRECTE)==-1)
	exit(-1);

if(fft(donnees)==-1)
	{
	ferme_fft();
	exit(-1);
	}

if(argc==4)
	{
	for(f=1;f>0;f-=0.05)
		{
		for(i=0;i<taille*2;i+=2)
			if(module(donnees[i],donnees[i+1],taille)>f)
				printf("*");
			else
				printf(" ");
		printf("\n");
		}
	}
printf("SIGNAL EN FREQUENCE\n");	
ferme_fft();
free(donnees);
exit(0);
}

#endif
