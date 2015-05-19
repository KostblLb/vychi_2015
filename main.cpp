#include <iostream>
#include <string>
#include <cmath>

double delta;
int m; 
const double a = 0.8, b = 1.2;
double** U;
double** U_new;
double** U_new1;
double** U_norm;
double lambda_max = -1.0, lambda_min = 1.0, lambda_max_old = 0.0, lambda_min_old = 0.0;
double** A(double** u, double a, double b){
	double** u_new = new double*[m];
	for (int i = 0; i < m; i++)
		u_new[i] = new double[m];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++){
			u_new[i][j] = u[i][j];
		}
	}
	for (int i = 1; i < m - 1; i++){
		for (int j = 1; j < m - 1; j++){
			u_new[i][j] = -(a)*(u[i][j - 1] + u[i][j + 1] - (2.0*u[i][j])) - (b)*(u[i - 1][j] + u[i + 1][j] - (2.0*u[i][j]));

		}
	}
	return u_new;
}
double** Amod(double** u, double a, double b){
	double** u_new = new double*[m];
	for (int i = 0; i < m; i++)
		u_new[i] = new double[m];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++){
			u_new[i][j] = u[i][j];
		}
	}
	for (int i = 1; i < m / 2; i++){
		for (int j = 1; j < m - 1; j++){
			u_new[i][j] = -(a)*(u[i][j - 1] + u[i][j + 1] - (2.0*u[i][j])) - (b)*(u[i - 1][j] + u[i + 1][j] - (2.0*u[i][j]));
			u_new[i][j] *= m*m;
		}
	}
	for (int j = m / 2 + 1; j < m - 1; j++){
		u_new[m / 2][j] = -(a)*(u[m / 2][j - 1] + u[m / 2][j + 1] - (2.0*u[m / 2][j])) - (b)*(u[m / 2 - 1][j] + u[m / 2 + 1][j] - (2.0*u[m / 2][j]));
		u_new[m / 2][j] *= m*m;
	}
	for (int i = m / 2 + 1; i < m - 1; i++){
		for (int j = m / 2 + 1; j < m - 1; j++)
			if (j > m / 2 && i < m*1.5 - j - 1){
				u_new[i][j] = -(a)*(u[i][j - 1] + u[i][j + 1] - (2.0*u[i][j])) - (b)*(u[i - 1][j] + u[i + 1][j] - (2.0*u[i][j]));
				u_new[i][j] *= m*m;
			}
	}
	return u_new;
}

double** B(double** u, double a, double b){
	double** u_new = new double*[m];
	for (int i = 0; i < m; i++)
		u_new[i] = new double[m];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++){
			u_new[i][j] = u[i][j];
		}
	}
	//double** u_new = A(u, a, b);
	for (int i = 1; i < m - 1; i++){
		for (int j = 1; j < m - 1; j++){
			u_new[i][j] = lambda_max*u[i][j] / (m*m) - (-(a)*(u[i][j - 1] + u[i][j + 1] - (2.0*u[i][j])) - (b)*(u[i - 1][j] + u[i + 1][j] - (2.0*u[i][j])));
		}
	}
	return u_new;
}
double mult(double** a, double** b){
	double s = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++){
			s += a[i][j] * b[i][j];
		}
	return s;
}
double** mult(double** a, double b){
	double** res = new double*[m];
	for (int i = 0; i < m; i++){
		res[i] = new double[m];
	}
	for (int i = 0; i < m; i++){
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j] * b;
	}
	return res;
}

double phi(double x, double y){
	return (x*x / (m*m))*sin(y / m);
}
double f(double x, double y){
	return 0.4*sin(y / m)*((3.0 * x*x / (m*m)) - 4.0);
}
double norma(double* x){
	double n = 0;
	for (int i = 0; i < m; i++){
		n += x[i] * x[i];
	}
	return sqrt(n);
}
double norma(double** x){
	return sqrt(mult(x, x));
}

double** subs(double** a, double** b){
	double** res = new double*[m];
	for (int i = 0; i < m; i++){
		res[i] = new double[m];
	}
	for (int i = 0; i < m; i++){
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j] - b[i][j];
	}
	return res;
}
double** subs(double** a, double b(double, double)){
	double** res = new double*[m];
	for (int i = 0; i < m; i++){
		res[i] = new double[m];
	}
	for (int i = 0; i < m; i++){
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j] - b(j, i);
	}
	return res;
}
double** subsmod(double** a, double** b){
	double** res = new double*[m];
	for (int i = 0; i < m; i++){
		res[i] = new double[m];
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j];
	}
	for (int i = 1; i < m / 2; i++){
		for (int j = 1; j < m - 1; j++){
			res[i][j] = a[i][j] - b[i][j];
		}
	}
	for (int j = m / 2 + 1; j < m - 1; j++)
		res[m / 2][j] = a[m / 2][j] - b[m / 2][j];
	for (int i = m / 2 + 1; i < m - 1; i++){
		for (int j = m / 2 + 1; j < m - 1; j++)
			if (j > m / 2 && i < m*1.5 - j - 1)
				res[i][j] = a[i][j] - b[i][j];
	}
	return res;
}
double** subsmod(double** a, double b(double, double)){
	double** res = new double*[m];
	for (int i = 0; i < m; i++){
		res[i] = new double[m];
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j];
	}
	for (int i = 1; i < m / 2; i++){
		for (int j = 1; j < m - 1; j++){
			res[i][j] = a[i][j] - b(j, i);
		}
	}
	for (int j = m / 2 + 1; j < m - 1; j++)
		res[m / 2][j] = a[m / 2][j] - b(j, m / 2);
	for (int i = m / 2 + 1; i < m - 1; i++){
		for (int j = m / 2 + 1; j < m - 1; j++)
			if (j > m / 2 && i < m*1.5 - j - 1)
				res[i][j] = a[i][j] - b(j, i);
	}
	return res;
}
void deletearr(double** arr){
	if (!arr) return;
	for (int i = 0; i < m; i++){
		delete[] arr[i];
	}
	delete[] arr;
}
int main(){
	std::cout << "enter table width and delta\n";
	std::cin >> m >> delta;
	U = new double*[m];
	//U_norm = new double*[m]; //
	for (int i = 0; i < m; i++){
		U[i] = new double[m];
		//U_norm[i] = new double[m]; //
	}
<<<<<<< HEAD
	for (int i = 0; i < m; i++){ //заполняем весь массив единицами
=======
	for (int i = 0; i < m; i++){ //çàïîëíÿåì âåñü ìàññèâ åäèíèöàìè
>>>>>>> origin/master
		for (int j = 0; j < m; j++){
			U[i][j] = 1;
		}
	}
	for (int i = 0; i < m; i++){
		U[i][0] = 0;
		U[0][i] = 0;
		U[i][m - 1] = 0;
		U[m - 1][i] = 0;
	}
<<<<<<< HEAD
	for (int i = m / 2; i < m; i++){//заполняем выбранные части массива нулями
=======
	for (int i = m / 2; i < m; i++){//çàïîëíÿåì âûáðàííûå ÷àñòè ìàññèâà íóëÿìè
>>>>>>> origin/master
		for (int j = 0; j < m; j++){
			if (j <= m / 2 || i >= (1.5*m) - j)
				U[i][j] = 0;
		}
	}
	double delta1 = 1e-1;
	double delta2 = 1e-4;
	double delta3 = 1e-5;
	double delta4 = 1e-10;
	int n = 0;
	std::string cause;
	while (true){

		lambda_max_old = lambda_max;
		U_new = A(U, a, b);
		lambda_max = mult(U_new, U);
		for (int i = 0; i < m; i++){
			delete[] U[i];
		}
		delete[] U;
		U = U_new;
		double norma = sqrt(mult(U, U)) / m;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				U[i][j] /= norma;
		n++;
		if (abs((lambda_max - lambda_max_old) / lambda_max) < delta2) {
<<<<<<< HEAD
			cause = "Converged";  //случай сходимости
			break;
		}
		if (lambda_max == 0.0) {
			cause = "NaN"; lambda_max = lambda_max_old;  //случай переполнения double
=======
			cause = "Converged";  //ñëó÷àé ñõîäèìîñòè
			break;
		}
		if (lambda_max == 0.0) {
			cause = "NaN"; lambda_max = lambda_max_old;  //ñëó÷àé ïåðåïîëíåíèÿ double
>>>>>>> origin/master
			break;
		}
	}
	for (int i = 0; i < m; i++){
		delete[] U[i];
		//		delete[] U_norm[i];
	}
	delete[] U;
	//	delete[] U_norm;
	std::cout << cause << ' ' << lambda_max << ' ' << n << '\n';
	system("pause");


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
<<<<<<< HEAD
	U = new double*[m]; ///всё сбрасываем для поиска lambda_min
=======
	U = new double*[m]; ///âñ¸ ñáðàñûâàåì äëÿ ïîèñêà lambda_min
>>>>>>> origin/master
	//U_norm = new double*[m]; //
	for (int i = 0; i < m; i++){
		U[i] = new double[m];
		//U_norm[i] = new double[m]; //
	}
<<<<<<< HEAD
	for (int i = 0; i < m; i++){ //заполняем весь массив единицами
=======
	for (int i = 0; i < m; i++){ //çàïîëíÿåì âåñü ìàññèâ åäèíèöàìè
>>>>>>> origin/master
		for (int j = 0; j < m; j++){
			U[i][j] = 1;
		}
	}
	for (int i = 0; i < m; i++){
		U[i][0] = 0;
		U[0][i] = 0;
		U[i][m - 1] = 0;
		U[m - 1][i] = 0;
	}
<<<<<<< HEAD
	for (int i = m / 2; i < m; i++){//заполняем выбранные части массива нулями
=======
	for (int i = m / 2; i < m; i++){//çàïîëíÿåì âûáðàííûå ÷àñòè ìàññèâà íóëÿìè
>>>>>>> origin/master
		for (int j = 0; j < m; j++){
			if (j <= m / 2 || i >= (1.5*m) - j)
				U[i][j] = 0;
		}
	}
	n = 0;
	while (true){
		lambda_min_old = lambda_min;
		double norma = sqrt(mult(U, U)) / m;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				U[i][j] /= norma;
		//norm(U);
		//U_new1 = B(U_norm, a, b);
		//U_new1 = B(U, a, b);
		U_new = B(U, a, b);
		lambda_min = mult(U_new, U);
		for (int i = 0; i < m; i++){
			delete[] U[i];
			//delete[] U_new1[i];
		}
		delete[] U;
		//delete[] U_new1;
		U = U_new;
		n++;
		if (abs((lambda_min - lambda_min_old) / lambda_min) < delta2) {
			cause = "Converged";
			break;
		}
		if (lambda_min == 0.0) {
			cause = "NaN"; lambda_min = lambda_min_old;
			break;
		}
	}
	for (int i = 0; i < m; i++){
		delete[] U[i];
		//		delete[] U_norm[i];
	}
	delete[] U;
	//	delete[] U_norm;
	lambda_min = (lambda_max - lambda_min);
	std::cout << cause << ' ' << lambda_min << ' ' << n << '\n';
	system("pause");


<<<<<<< HEAD
	//поиск оператора по методу минимальных невязок

	U = new double*[m]; ///всё сбрасываем для поиска lambda_min
=======
	//ïîèñê îïåðàòîðà ïî ìåòîäó ìèíèìàëüíûõ íåâÿçîê

	U = new double*[m]; ///âñ¸ ñáðàñûâàåì äëÿ ïîèñêà lambda_min
>>>>>>> origin/master
	//U_norm = new double*[m]; //
	for (int i = 0; i < m; i++){
		U[i] = new double[m];
		//U_norm[i] = new double[m]; //
	}
<<<<<<< HEAD
	for (int i = 0; i < m; i++){ //заполняем весь массив нулями(для нач. приближения)
=======
	for (int i = 0; i < m; i++){ //çàïîëíÿåì âåñü ìàññèâ íóëÿìè(äëÿ íà÷. ïðèáëèæåíèÿ)
>>>>>>> origin/master
		for (int j = 0; j < m; j++){
			U[i][j] = 0;
		}
	}
	for (int i = 0; i < m / 2; i++)
		U[i][m - 1] = phi(m - 1, i);
<<<<<<< HEAD
	for (int i = m / 2; i < m; i++){//заполняем границу по граничным условиям
=======
	for (int i = m / 2; i < m; i++){//çàïîëíÿåì ãðàíèöó ïî ãðàíè÷íûì óñëîâèÿì
>>>>>>> origin/master
		for (int j = 0; j < m; j++){

			if (j < m / 2 || i >(1.5*m) - j)
				U[i][j] = 0;
			if (j == m / 2 || i == m*1.5 - j - 1 || (i == m / 2 && j < m / 2))
				U[i][j] = phi(j, i);
		}
	}
	double tau = 0, **r = new double*[m], **z = NULL;
	for (int i = 0; i < m; i++){
		r[i] = new double[m];
		for (int j = 0; j < m; j++){
			r[i][j] = 1;
		}
	}
	n = 0;
	double ** Uex = new double*[m];
	for (int i = 0; i < m; i++) Uex[i] = new double[m];
<<<<<<< HEAD
	for (int i = 0; i < m; i++){ //заполняем весь массив нулями(для нач. приближения)
=======
	for (int i = 0; i < m; i++){ //çàïîëíÿåì âåñü ìàññèâ íóëÿìè(äëÿ íà÷. ïðèáëèæåíèÿ)
>>>>>>> origin/master
		for (int j = 0; j < m; j++){
			Uex[i][j] = 0;
		}
	}
	for (int i = 0; i < m / 2; i++)
		for (int j = 0; j < m; j++)
			Uex[i][j] = phi(j, i);
<<<<<<< HEAD
	for (int i = m / 2; i < m; i++){//заполняем границу по граничным условиям
=======
	for (int i = m / 2; i < m; i++){//çàïîëíÿåì ãðàíèöó ïî ãðàíè÷íûì óñëîâèÿì
>>>>>>> origin/master
		for (int j = 0; j < m; j++){

			if (j < m / 2 || i >(1.5*m) - j)
				Uex[i][j] = 0;
			if ((j >= m / 2 && i <= m*1.5 - j - 1) || (i == m / 2 && j < m / 2))
				Uex[i][j] = phi(j, i);
		}
	}
	double normUex = norma(Uex);
	double normr1 = 1, normr = 0, normz, normU = 0, normU1 = 0;
	double nr = 1, nr1 = 0;

	while (true){
		//1)rk=Ayk-f
		double** Ayk = Amod(U, a, b);
		for (int i = 0; i < m / 2; i++)
			//r[i][m - 1] = Ayk[i][m-1]-f(m-1,i);
			Ayk[i][m - 1] = 0;
<<<<<<< HEAD
		for (int i = m / 2; i < m; i++){//заполняем границу по граничным условиям
=======
		for (int i = m / 2; i < m; i++){//çàïîëíÿåì ãðàíèöó ïî ãðàíè÷íûì óñëîâèÿì
>>>>>>> origin/master
			for (int j = 0; j < m; j++){
				if (j < m / 2 || i >(1.5*m) - j)
					Ayk[i][j] = 0;
				if (j == m / 2 || i == m*1.5 - j - 1 || (i == m / 2 && j < m / 2))
					//r[i][j] = Ayk[i][j] - f(j,i);
					Ayk[i][j] = 0;
			}
		}
		deletearr(r);
		r = subsmod(Ayk, f);

		for (int i = 0; i < m / 2; i++)
			//r[i][m - 1] = Ayk[i][m-1]-f(m-1,i);
			r[i][m - 1] = 0;
<<<<<<< HEAD
		for (int i = m / 2; i < m; i++){//заполняем границу по граничным условиям
=======
		for (int i = m / 2; i < m; i++){//çàïîëíÿåì ãðàíèöó ïî ãðàíè÷íûì óñëîâèÿì
>>>>>>> origin/master
			for (int j = 0; j < m; j++){
				if (j < m / 2 || i >(1.5*m) - j)
					r[i][j] = 0;
				if (j == m / 2 || i == m*1.5 - j - 1 || (i == m / 2 && j < m / 2))
					//r[i][j] = Ayk[i][j] - f(j,i);
					r[i][j] = 0;
			}
		}
		deletearr(Ayk);
		//2)tau=(Ar,r)/(Ar,Ar)
		double** Ar = Amod(r, a, b);
		tau = mult(Ar, r) / mult(Ar, Ar);
		//3)yk+1=y-tau*r
		double** tr = mult(r, tau);
		U_new = subsmod(U, tr);
		deletearr(Ar);
		deletearr(tr);
		normU1 = normU;
		normU = norma(U_new);
		
		deletearr(U);
		U = U_new;
		n++;
		if (abs(normU-normU1) < delta){
			z = subsmod(U, phi);
			for (int i = 0; i < m / 2; i++)
				z[i][m - 1] = 0;
			for (int i = m / 2; i < m; i++){
				for (int j = 0; j < m; j++){
					if (j < m / 2 || i >(1.5*m) - j)
						z[i][j] = 0;
					if (j == m / 2 || i == m*1.5 - j - 1 || (i == m / 2 && j < m / 2))
						z[i][j] = 0;
				}
			}
			normz = norma(z);
			nr = norma(r);
			fprintf(stdout, "w:%d  delta:%f  n:%d  |U|:%f  |z|:%f  |r|:%f\n", m, delta, n, normU, normz, nr);
<<<<<<< HEAD


=======
>>>>>>> origin/master
			for (int re = m-1; re >=0; --re) {
				for (int rere = 0; rere < m; ++rere)
					fprintf(stdout, "%f ", U[re][rere]);
				fprintf(stdout, "\n");
			}
			fprintf(stdout, "--------------------------\n");

			for (int re = m-1; re >=0; --re) {
				for (int rere = 0; rere < m; ++rere)
					fprintf(stdout, "%f ", Uex[re][rere]);
				fprintf(stdout, "\n");
			}
<<<<<<< HEAD

=======
>>>>>>> origin/master
			break;
		}
	}

	deletearr(r);
	deletearr(U);
	deletearr(Uex);
	system("pause");
}
