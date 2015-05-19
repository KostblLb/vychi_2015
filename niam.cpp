//#include <iostream>
//#include <string>
//int m;
//const double a = 0.8, b = 1.2;
//double delta;
//
//void deletearr(double** arr){
//	if (!arr) return;
//	for (int i = 0; i < m; i++){
//		delete[] arr[i];
//	}
//	delete[] arr;
//}
//
//void fill(double**& mat, double num){
//	deletearr(mat);
//	mat = new double*[m];
//	for (int i = 0; i < m; i++){
//		mat[i] = new double[m];
//		for (int j = 0; j < m; j++)
//			mat[i][j] = num;
//	}
//}
//void fill(double**& mat, double func(double, double)){
//	deletearr(mat);
//	mat = new double*[m];
//	for (int i = 0; i < m; i++){
//		mat[i] = new double[m];
//		for (int j = 0; j < m; j++)
//			mat[i][j] = func(j, i);
//	}
//}
//
//void fill(double**& mat, double** other){
//	deletearr(mat);
//	mat = new double*[m];
//	for (int i = 0; i < m; i++){
//		mat[i] = new double[m];
//		for (int j = 0; j < m; j++)
//			mat[i][j] = other[i][j];
//	}
//}
//
//double** A(double** u){
//	double** u_new = new double*[m];
//	for (int i = 0; i < m; i++)
//		u_new[i] = new double[m];
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < m; j++){
//			u_new[i][j] = u[i][j];
//		}
//	}
//	for (int i = 1; i < m - 1; i++){
//		for (int j = 1; j < m - 1; j++){
//			u_new[i][j] = -(a)*(u[i][j - 1] + u[i][j + 1] - (2.0*u[i][j])) - (b)*(u[i - 1][j] + u[i + 1][j] - (2.0*u[i][j]));
//			u_new[i][j] *= m*m;
//		}
//	}
//	return u_new;
//}
//
//double mult(double** a, double** b){
//	double s = 0;
//	for (int i = 0; i < m; i++)
//		for (int j = 0 ; j < m; j++){
//			s += a[i][j] * b[i][j];
//		}
//	return s;
//}
//
//double** mult(double** a, double b){
//	double** res = NULL;
//	fill(res, a);
//	for (int i = 0; i < m; i++){
//		for (int j = 0; j < m; j++)
//			res[i][j] *= b;
//	}
//	return res;
//}
//
//double phi(double x, double y){
//	return (x*x / (m*m))*sin(y / m);
//}
//
//double f(double x, double y){
//	return 0.4*sin(y / m)*((3.0 * x*x / (m*m)) - 4.0);
//}
//
//double norma(double* x){
//	double n = 0;
//	for (int i = 0; i < m; i++){
//		n += x[i] * x[i];
//	}
//	return sqrt(n);
//}
//
//double norma(double** x){
//	return sqrt(mult(x, x));
//}
//
//double** subs(double** a, double** b){
//	double** res = NULL;
//	fill(res, a);
//	for (int i = 1; i < m-1; i++){
//		for (int j = 1; j < m-1; j++)
//			res[i][j]  -= b[i][j];
//	}
//	return res;
//}
//
//double** subs(double** a, double b(double, double)){
//	double** res = NULL;
//	fill(res, a);
//	for (int i = 1; i < m-1; i++){
//		for (int j = 1; j < m-1; j++)
//			res[i][j] -= b(j, i);
//	}
//	return res;
//}
//
//void display(double** mat){
//	for (int i = m - 1; i >= 0; i--){
//		for (int j = 0; j < m; j++){
//			std::cout << mat[i][j] << ' ';
//		}
//		std::cout << std::endl;
//	}
//}
//int main(){
//	std::cin >> m >> delta;
//	std::cout << std::endl;
//
//	double** U = NULL; //сетка для итерирования
//	fill(U, 0.0);
//	for (int i = 0; i < m; i++){
//		U[i][0] = phi(0, i);
//		U[0][i] = phi(i, 0);
//		U[i][m - 1] = phi(m - 1, i);
//		U[m - 1][i] = phi(i, m - 1);
//	}
//	double** Uex = NULL; //сетка точных значений
//	fill(Uex, phi);
//	double nuex = norma(Uex);
//	double nr1 = 1, nr2 = 0;
//	int n = 0;
//	FILE* file;
//	fopen_s(&file, "result.txt", "wb");
//	while (true){
//		double** AU = A(U);
//		double** r = subs(AU, f);
//		nr2 = nr1;
//		nr1 = norma(r);
////		display(r);
//		double** Ar = A(r);
//		double tau = mult(Ar, r) / mult(Ar, Ar);
//		double** tau_r = mult(r, tau);
//		double** Unew = subs(U, tau_r);
////		display(Unew);
//		double nu1 = norma(U);
//		double nu2 = norma(Unew);
//		double** z = subs(U, Uex);
//		double nz = norma(z);
//		deletearr(z);
//		deletearr(r);
//		deletearr(AU);
//		deletearr(Ar);
//		deletearr(tau_r);
//		deletearr(U);
//		U = Unew;
//		n++;
//		if (abs((nr1-nr2) / nr2) < delta){
//			std::cout << "\n\n";
//			display(U);
//			std::cout << "\n\n";
//			display(Uex);
//			std::cout << n << '\n';
//			fprintf(file, "w:%d  delta:%f  n:%d  |U|:%f  |z|:%f  |r|:%f\n", m, delta, n, nu2, nz, nr1);
//			fprintf(stdout, "w:%d  delta:%f  n:%d  |U|:%f  |z|:%f  |r|:%f\n", m, delta, n, nu2, nz, nr1);
//			break;
//		}
//	}
//	fclose(file);
//	system("pause");
//}