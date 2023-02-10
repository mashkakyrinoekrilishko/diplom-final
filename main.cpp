#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <set>
#include <algorithm>
#include <iterator>

using namespace std;

struct node //структура узла для 2D сетки
{
	double r;
	double z;
};

struct node_3D //структура узла для 3D сетки
{
	double x;
	double y;
	double z;
};

struct edge // структура ребра
{
	node_3D beg, end, middle;
};

struct material
{
	double mu;
	double lambda;
	double sigma;
};

struct element
{
	vector<int> num; //локальные узлы
	int mater;
};

struct layer {
	double y0, y1;
	unsigned int n_mat;
};

static const double G1[4][4] =
{
	{2.0, 1.0, -2.0, -1.0},
	{1.0, 2.0, -1.0, -2.0},
	{-2.0, -1.0, 2.0, 1.0},
	{-1.0, -2.0, 1.0, 2.0}
};

static const double G2[4][4] =
{
	{2.0, -2.0, 1.0, -1.0},
	{-2.0, 2.0, -1.0, 1.0},
	{1.0, -1.0, 2.0, -2.0},
	{-1.0, 1.0, -2.0, 2.0}
};

static const double G3[4][4] =
{
	{-2.0, 2.0, -1.0, 1.0},
	{-1.0, 1.0, -2.0, 2.0},
	{2.0, -2.0, 1.0, -1.0},
	{1.0, -1.0, 2.0, -2.0}
};

static const double G3T[4][4] =
{
	{-2.0, -1.0, 2.0, 1.0},
	{2.0, 1.0, -2.0, -1.0},
	{-1.0, -2.0, 1.0, 2.0},
	{1.0, 2.0, -1.0, -2.0}
};

static const double M[4][4] =
{
	{4.0, 2.0, 2.0, 1.0},
	{2.0, 4.0, 1.0, 2.0},
	{2.0, 1.0, 4.0, 2.0},
	{1.0, 2.0, 2.0, 4.0}
};


vector<node> nodes;         // все узлы в порядке глобальной нумерации
vector<element> elems;      // все элементы в прядке глобальной нумерации
vector<material> materials; // все материалы по индексам
vector<int>KR1;    // KR1[i][j] на j-ом узле заданы i-ые краевые 1 рода
vector<int> KR1_3D;    // KR1[i][j] на j-ом узле заданы i-ые краевые 1 рода
vector<layer> layers;       //слои в среде
vector<edge> edges;         // все ребра в порядке глобальной нумерации
vector<element> obj_elems;      // все элементы объекта в прядке глобальной нумерации

class KURS
{
public:
	vector<double> di, al, au, b_loc, b, time, M2, M1, temp, p, z, r, Ar, y, L, D, U, ar, az; //x0, x1, x2, 
	vector<int> ja, ia;
	vector<node> Arz;
	vector<vector<double>> A_loc, G_loc, M_loc, M_loc_g, x0, x0_3D;
	set <double> r_s, z_s; //границы КЭ
	int N, Kel, time_number, istoc, num_r, num_z, n_layers, maxIter = 100000, iter = 0, Nmat, rasm; //rasm = 4 или 12 в зависимости от 2 и 3 дименшн
	double eps = 1E-6, eps_time = 1E-16, normR = 0, normB, hr, hz, rp, lambda, sigma, zs, t0, t1, t2, nu0, nu1, nu2, istoc_r, istoc_z, h_r, h_z, d_r, d_z, r_left, r_right, z_top, z_bottom, usel;

	///----------------общие функции и процедуры--------------------------------------------
	double gamma(double r); //гамма 
	void Input(); //считывание всех необходимых данных
	void Clear_and_resize(); //инициализация векторов
	void multiply(vector<vector<double>>& A, vector<double>& x, vector<double>& res);

	//------------------для двумерной задачи--------------------------------------------------
	void Generate_Portrait(); // генерация портрета
	void Assemble_Locals(int el_id); // внесение локальных в глобальную СЛАУ
	void Generate_2D_Net(string net_file); //генерация сетки и конечных элементов для 2-мерной задачи
	void Read_Layers(string obj_file); //формируем слоистую среду
	void Get_G(); //получение локальной матрицы G
	void Get_M();//получение локальных матриц M
	void Locals(int el_id, double t); // получение локальной матрицы А
	void Find_field(); //цикл по времени
	void Get_KR1(double t); // учет первых краевых
	int  Get_Num_Layer(double x0, double y0, double x1, double y1); //определяем к какому слою относится элемент

	//---------------------для трехмерной задачи-----------------------------------------------
	double Get_solution(int time_layer,double r, double z, double variable, bool isB);
	void Generate_3D_Net(string net_file, bool go_to_z_up); //генерация сетки и конечных элементов
	void Middles();
	void Find_in_object(); //цикл по времени
	void Locals_3D(int el_id, int time_layer);//dsx
	void Get_G_3D(); //получение локальной матрицы G
	void Get_M_3D();//получение локальных матриц M
	void Get_b_3D(int time_layer, int el_id); // получение локального b
	void Assemble_Locals_3D(int el_id); // внесение локальных A, b  в глобальную СЛАУ
	void Generate_Portrait_3D();
	int Get_Num_Layer_3D(double z0, double z1);
	double Get_solution_3D(double t, double x, double y, double z, int i, bool B);
	void Make_grid(bool go_to_z_up);
	void Make_edges();
	void Make_obj_elems();
	void Get_KR1_3D();
	void Make_centers();
	void diag_preconditioning(vector<double>& x0);
	double h_x, h_y, h_obj_z, d_x, d_y, d_obj_z, x_left, x_right, y_left, y_right, z_obj_top, z_obj_bottom, hx, hy, 
		hz_o, x_obj_left, x_obj_right, h_obj_x, d_obj_x, y_obj_left, y_obj_right,  h_obj_y, d_obj_y, z_sr_left, z_sr_right, d_sr_z, reciever_x, reciever_y, reciever_z, time_h, time_beg, time_d;
	int num_x, num_y, num_obj_z, num_edges, num_elems, lay_obj;
	vector<double> ax, ay, a_obj_z, F, otvet;
	vector <node_3D> centers;
	set <double> x_s, y_s, z_obj_s; //границы КЭ

	//---------------------для решателя--------------------------------------------------------
	void LOS_LU(int t);
	void LOS_LU_3D(int t);
	void FactLU(vector<double>& L, vector<double>& U, vector<double>& D);
	void Direct(vector<double>& L, vector<double>& D, vector<double>& y, vector<double>& b);
	void Reverse(vector<double>& U, vector<double>& x, vector<double>& y);
	double Norm(vector<double>& x);
	double mult(const vector<double>& a, const vector<double>& b);
	void Ax(vector<double>& x, vector<double>& y);
	void ATx(vector<double>& x, vector<double>& y);
	double relative_residual(vector<double>& x);
	void llt_preconditioning(vector<double>& x0);
	void llt();
	void solve_auxiliary_system(const vector<double>& f, vector<double>& x);
	void vec_diff(const vector<double>& x, const vector<double>& y, vector<double>& res);
};

//--------------------блок общих функций и процедур---------------------
double KURS::gamma(double r) // значение гамма всегда равно 1/r^2
{
	return 1.0 / (r * r);
}

void KURS::Input() // чтение данных
{
	ifstream in;

	in.open("istoc.txt");
	in >> istoc_r >> istoc_z;
	in.close();

	Read_Layers("layers");
	Generate_2D_Net("gen_net");

	for (int i = 0; i < N; i++)
		if (nodes[i].r == r_right || nodes[i].z == z_top || nodes[i].z == z_bottom)
			KR1.push_back(i);

	in.open("material.txt");
	in >> Nmat;
	materials.resize(Nmat + 1);
	for (int i = 0; i < Nmat; i++)
	{
		in >> materials[i].mu;
		materials[i].lambda = 1. / materials[i].mu;
		in >> materials[i].sigma;
	}

	in.close();

	in.open("obj_material.txt"); 	//задаем параметры материалов для объекта 
	in >> materials[Nmat].mu >> materials[Nmat].sigma;
	materials[Nmat].lambda = 1. / materials[Nmat].mu;
	in.close();

	in.open("time.txt");
	in >> time_number >> time_beg >> time_h >> time_d; //считываем количество временных слоев, начало по времени, шаг и коэф
	time.resize(time_number);
	time[0] = time_beg;
	/*time_d = sqrt(time_d);
	time_h = time_h / (1 + time_d);*/
	for (int i = 1; i < time_number; i++) //считываем временные слои
		time[i] = time[i-1] +time_h*pow(time_d,i-1);
	in.close();

	Generate_3D_Net("gen_obj", false); //генерируем сетку для объекта
	Middles();

	in.open("reciever.txt");
	in >> reciever_x >> reciever_y >> reciever_z;
	in.close();
}

void KURS::Clear_and_resize()
{
	au.clear();
	au.resize(ia[N]);
	al.clear();
	al.resize(ia[N]);
	di.clear();
	di.resize(N);
	b.clear();
	b.resize(N);
	ja.resize(ia[N]);
	G_loc.resize(rasm);
	M_loc.resize(rasm);
	A_loc.resize(rasm);
	b_loc.clear();
	b_loc.resize(rasm);
	M_loc_g.resize(rasm);
	for (int i = 0; i < rasm; i++)
	{
		A_loc[i].resize(rasm);
		M_loc[i].resize(rasm);
		M_loc_g[i].resize(rasm);
		G_loc[i].resize(rasm);
	}
	M1.resize(N);
	M2.resize(N);
}

void KURS::multiply(vector<vector<double>>& A, vector<double>& x, vector<double>& res)
{
	for (int i = 0; i < rasm; i++)
	{
		res[i] = 0.;
		for (int j = 0; j < rasm; j++)
			res[i] += A[i][j] * x[j];
	}
}
//----------------------------------------------------------------------


//-------------блок функций и процедур для 2D задачи--------------------
void KURS::Generate_2D_Net(string net_file)
{
	ifstream fin(net_file + ".txt");
	fin >> r_left >> r_right >> z_top >> z_bottom >> h_r >> d_r >> h_z >> d_z; //считываем границы области, шаги по r и z, коэффициент растяжения

   //строим сетку по r
	num_r = 0; //обнуляем число элементов по r 
	double h = h_r / d_r;
	bool flag = true;
	ar.push_back(istoc_r);
	double usel = istoc_r;
	for (int i = 0; usel > r_left; i++) //идем влево от источника
	{
		flag = true;
		h = h * d_r;
		usel = ar[i] - h;
		if (usel < r_left)
		{
			flag = false;
			usel = r_left;
			if (abs(ar[i] - usel) < abs(ar[i] - ar[i - 1]))
			{
				ar[i] = usel;
				h = abs(ar[i] - ar[i - 1]) / d_r;
				i--;
			}
			else
				ar.push_back(usel);
		}
		if (flag == true)
			ar.push_back(usel);
	}

	reverse(ar.begin(), ar.end());

	h = h_r / d_r;
	for (int i = ar.size() - 1; usel < r_right; i++) //идем вправо от источника
	{
		flag = true;
		h = h * d_r;
		usel = ar[i] + h;
		if (usel > r_right)
		{
			usel = r_right;
			flag = false;
			if (abs(ar[i] - usel) < abs(ar[i] - ar[i - 1]))
			{
				ar[i] = usel;
				h = abs(ar[i] - ar[i - 1]) / d_r;
				i--;
			}
			else
				ar.push_back(usel);
		}
		if (flag == true)
			ar.push_back(usel);
	}

	num_r = ar.size();

	for (int j = 0; j < num_r; j++)
		r_s.insert(ar[j]);

	//строим сетку по z 
	num_z = 0; //обнуляем число элементов по z

	az.push_back(istoc_z);
	usel = istoc_z;
	h = h_z / d_z;
	for (int i = 0; usel > z_bottom; i++) //идем вниз от источника
	{
		flag = true;
		h = h * d_z;
		usel = az[i] - h;
		if (i != 0)
		{
			for (int j = 0; j < n_layers; j++) //проеряем значения в окрестностях границ слоев
				if (layers[j].y0 > usel && az[i] > layers[j].y0)
				{
					flag = false;
					usel = layers[j].y0;
					if (abs(az[i] - usel) < abs(az[i] - az[i - 1]))
					{
						az[i] = usel;
						h = abs(az[i] - az[i - 1]) / d_z;
						i--;
					}
					else
						az.push_back(usel);
				}
		}
		if (flag == true)
			az.push_back(usel);
	}

	reverse(az.begin(), az.end());
	h = h_z / d_z;
	for (int i = az.size() - 1; usel < z_top; i++) //идем вверх от источника
	{
		flag = true;
		h = h * d_z;
		usel = az[i] + h;
		if (usel > z_top)
		{
			usel = z_top;
			flag = false;
			if (abs(az[i] - usel) < abs(az[i] - az[i - 1]))
			{
				az[i] = usel;
				h = abs(az[i] - az[i - 1]) / d_z;
				i--;
			}
			else
				az.push_back(usel);
		}
		if (flag == true)
			az.push_back(usel);
	}

	num_z = az.size();

	//обработка по z
	for (int j = 0; j < num_z; j++)
		z_s.insert(az[j]);

	//для вывода разбиений по r и z для построения сетки в питоне
    for (int i = 0; i < num_r; i++)
		cout << ar[i] << ", ";
	cout << endl;

	for (int i = 0; i < num_z; i++)
		cout << az[i] << ", ";

		// Формирование узлов
	Kel = (num_r - 1) * (num_z - 1); //количество конечных элементов 
	N = num_r * num_z; //число узлов
	nodes.resize(N);
	int k = 0;
	for (set<double>::iterator j = z_s.begin(); j != z_s.end(); j++)
		for (set<double>::iterator i = r_s.begin(); i != r_s.end(); i++)
		{
			nodes[k].r = *i;
			nodes[k].z = *j;
			if (nodes[k].r == istoc_r && nodes[k].z == istoc_z)
				istoc = k;
			k++;
		}

	// Получаем КЭ
	int ku = 0;
	elems.resize(Kel);
	for (int i = 0; i < Kel; i++)
		elems[i].num.resize(4);
	set<double>::iterator j_it = z_s.begin();
	for (unsigned int j = 0; j < num_z - 1; j++)
	{
		set<double>::iterator i_it = r_s.begin();
		for (unsigned int i = 0; i < num_r - 1; i++)
		{
			elems[ku].num[0] = num_r * j + i;
			elems[ku].num[1] = num_r * j + i + 1;
			elems[ku].num[2] = num_r * (j + 1) + i;
			elems[ku].num[3] = num_r * (j + 1) + i + 1;
			double x0 = *i_it, y0 = *j_it;
			set<double>::iterator k = j_it;
			k++;
			i_it++;
			double x1 = *i_it, y1 = *k;
			unsigned int num_area = Get_Num_Layer(x0, y0, x1, y1);
			elems[ku].mater = layers[num_area].n_mat;
			ku++;
		}
		j_it++;
	}
}

int KURS::Get_Num_Layer(double x0, double y0, double x1, double y1)
{
	for (int i = 0; i < n_layers; i++)
		if (layers[i].y0 <= y0 && layers[i].y1 >= y1)
			return i;
	return 0;
}

void KURS::Read_Layers(string obj_file)
{
	ifstream fin(obj_file + ".txt");
	layer tmp;
	fin >> n_layers;
	for (int i = 0; i < n_layers; i++)
	{
		fin >> tmp.y0 >> tmp.y1 >> tmp.n_mat;
		layers.push_back(tmp);
	}
}

void KURS::Get_G() // получение локальной G
{
	double a1 = (lambda * hz * rp) / (6 * hr),
		a2 = (lambda * hz) / (12),
		a3 = (lambda * hr * rp) / (6 * hz),
		a4 = (lambda * hr * hr) / (12 * hz);
	G_loc[0][0] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
	G_loc[0][1] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
	G_loc[0][2] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
	G_loc[0][3] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;

	G_loc[1][0] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
	G_loc[1][1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
	G_loc[1][2] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
	G_loc[1][3] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;

	G_loc[2][0] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
	G_loc[2][1] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
	G_loc[2][2] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
	G_loc[2][3] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;

	G_loc[3][0] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
	G_loc[3][1] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;
	G_loc[3][2] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
	G_loc[3][3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
}

void KURS::Get_M() // получение локальной М
{
	double g1 = gamma(rp), g2 = gamma(rp + hr);
	M_loc_g[0][0] = M_loc_g[2][2] = hr * hz / 4 * (
		g1 * (rp / 3 + hr / 15) +
		g2 * (rp / 9 + 2 * hr / 45));
	M_loc_g[0][2] = M_loc_g[2][0] = hr * hz / 12 * (
		g1 * (rp / 2 + hr / 10) +
		g2 * (rp / 6 + hr / 15));
	M_loc_g[0][3] = M_loc_g[3][0] = hr * hz / 12 * (
		g1 * (rp / 6 + hr / 15) +
		g2 * (rp / 6 + hr / 10));
	M_loc_g[1][0] = M_loc_g[3][2] = M_loc_g[2][3] = M_loc_g[0][1] =
		hr * hz / 4 * (
			g1 * (rp / 9 + 2 * hr / 45) +
			g2 * (rp / 9 + hr / 15));
	M_loc_g[1][1] = M_loc_g[3][3] = hr * hz / 4 * (
		g1 * (rp / 9 + hr / 15) +
		g2 * (rp / 3 + 4 * hr / 15));
	M_loc_g[1][2] = M_loc_g[2][1] = hr * hz / 12 * (
		g1 * (rp / 6 + hr / 15) +
		g2 * (rp / 6 + hr / 10));
	M_loc_g[1][3] = M_loc_g[3][1] = hr * hz / 12 * (
		g1 * (rp / 6 + hr / 10) +
		g2 * (rp / 2 + 2 * hr / 5));

	M_loc[0][0] = M_loc[2][2] = sigma * hr * hz / 4 * (4 * rp / 9 + hr / 9);
	M_loc[0][2] = M_loc[2][0] = sigma * hr * hz / 12 * (2 * rp / 3 + hr / 6);
	M_loc[0][3] = M_loc[3][0] = sigma * hr * hz / 12 * (rp / 3 + hr / 6);
	M_loc[1][0] = M_loc[3][2] = M_loc[2][3] = M_loc[0][1] = sigma * hr * hz / 4 * (2 * rp / 9 + hr / 9);
	M_loc[1][1] = M_loc[3][3] = sigma * hr * hz / 4 * (4 * rp / 9 + hr / 3);
	M_loc[1][2] = M_loc[2][1] = sigma * hr * hz / 12 * (rp / 3 + hr / 6);
	M_loc[1][3] = M_loc[3][1] = sigma * hr * hz / 12 * (2 * rp / 3 + hr / 2);
}

void KURS::Locals(int el_id, double t) // получение локальной матрицы А
{
	element el = elems[el_id];
	hr = nodes[el.num[1]].r - nodes[el.num[0]].r;
	hz = nodes[el.num[2]].z - nodes[el.num[0]].z;
	rp = nodes[el.num[0]].r;
	zs = nodes[el.num[0]].z;
	lambda = materials[el.mater].lambda;
	sigma = materials[el.mater].sigma;
	Get_G();
	Get_M();
}

void KURS::Generate_Portrait() // генерация портрета
{
	ia.resize(N + 1);
	ja.resize(rasm * rasm * Kel);
	vector<int> t1(rasm * rasm * Kel), t2(rasm * rasm * Kel), beg(N);
	int s = 0;
	for (int i = 0; i < N; i++)
		beg[i] = 0;
	for (int el = 0; el < Kel; el++)
	{
		for (int i = 0; i < rasm; i++)
		{
			int k = elems[el].num[i];
			for (int j = i + 1; j < rasm; j++)
			{
				int ind1 = k;
				int ind2 = elems[el].num[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				int iaddr = beg[ind2];
				if (iaddr == 0)
				{
					s++;
					beg[ind2] = s;
					t1[s] = ind1;
					t2[s] = 0;
				}
				else
				{
					while (t1[iaddr] < ind1 && t2[iaddr] > 0)
						iaddr = t2[iaddr];
					if (t1[iaddr] > ind1)
					{
						s++;
						t1[s] = t1[iaddr];
						t2[s] = t2[iaddr];
						t1[iaddr] = ind1;
						t2[iaddr] = s;
					}
					else if (t1[iaddr] < ind1)
					{
						s++;
						t2[iaddr] = s;
						t1[s] = ind1;
						t2[s] = 0;
					}
				}
			}
		}
	}

	ia[0] = 0;
	for (int i = 0; i < N; i++)
	{
		ia[i + 1] = ia[i];
		int a = beg[i];
		while (a != 0)
		{
			ja[ia[i + 1]] = t1[a];
			ia[i + 1]++;
			a = t2[a];
		}
	}
}

void KURS::Assemble_Locals(int el_id) // внесение локальных A, b  в глобальную СЛАУ
{
	vector<int> L = elems[el_id].num;
	int k = elems[el_id].num.size(); // размерность локальной матрицы
	for (int i = 0; i < k; i++)
		di[L[i]] += A_loc[i][i];

	for (int i = 0; i < rasm; i++)
		for (int j = 0; j < i; j++)
			for (int k = ia[L[i]]; k < ia[L[i] + 1]; k++)
				if (ja[k] == L[j])
				{
					al[k] += A_loc[i][j];
					au[k] += A_loc[j][i];
					k++;
					break;
				}
}

void KURS::Find_field()
{
	x0.resize(time_number);
	for (int i = 0; i < time_number; i++)
		x0[i].resize(N);
	temp.resize(N);

	//-------Инициализация решения на 0-м слое как решение стационарной задачи----------------
	Clear_and_resize(); //очищаем вектора для каждого временного слоя

	for (int i = 0; i < Kel; i++)
	{
		Locals(i, time[0]);

		for (int i = 0; i < rasm; i++)
			for (int j = 0; j < rasm; j++)
				A_loc[i][j] = G_loc[i][j] + M_loc_g[i][j];

		Assemble_Locals(i); //сборка левой части

		for (int j = 0; j < rasm; j++)  //сборка вектора правой части
			b[elems[i].num[j]] += b_loc[j];
	}

	b[istoc] = b[istoc] + 1. / 2. / 3.14 / nodes[istoc].r; //задаем точечный источник в узле сетки

	Get_KR1(time[0]); //учет первых краевых

	LOS_LU(0);

	ofstream f("res.txt");
	f << "t = 0" << endl;
	for (int i = 0; i < N; i++)
		f << nodes[i].r << "\t" << nodes[i].z << "\t" << x0[0][i] << endl;
	f << endl;

	ofstream out("result.txt");
	out << "t = " << time[0] << endl;
	for (int i = 1; i <= N / num_r; i++)
	{
		for (int j = N - num_r * i; j < N - num_r * (i - 1); j++)
			out << x0[0][j] << "\t";
		out << endl;
	}
	out << endl;
	//------------------------------------------------------------------------------------------
	ofstream fout("check_Bfi.txt");
	//-----------------------Инициализация решения на 1-м слое----------------------------------
	for (int t = 1; t < time_number; t++)
	{
	//int t = 1;
		Clear_and_resize(); //очищаем вектора для каждого временного слоя

		//считаем коэффициенты по времени для двухслойной неявной схемы
		t0 = time[t] - time[t-1]; //t0
		nu0 = 1 / t0;

		for (int i = 0; i < Kel; i++)
		{
			Locals(i, time[t]);
			for (int i = 0; i < rasm; i++)
				for (int j = 0; j < rasm; j++)
					A_loc[i][j] = G_loc[i][j] + M_loc[i][j] * nu0 + M_loc_g[i][j];

			Assemble_Locals(i); //сборка левой части

			for (int j = 0; j < 4; j++) //сборка вектора правой части
				temp[j] = x0[t-1][elems[i].num[j]];
			multiply(M_loc, temp, M2);//M*x2
			for (int j = 0; j < 4; j++)
				M2[j] = M2[j] * nu0;
			for (int j = 0; j < 4; j++)
				b[elems[i].num[j]] += b_loc[j] + M2[j];
		}

		Get_KR1(time[t]);

		LOS_LU(t);

		f << "t = " << time[t] << endl;
		for (int i = 0; i < N; i++)
			f << nodes[i].r << "\t" << nodes[i].z << "\t" << x0[t][i] << endl;
		f << endl;

		out << "t = " << time[t] << endl;
		for (int i = 1; i <= N / num_r; i++)
		{
			for (int j = N - num_r * i; j < N - num_r * (i - 1); j++)
				out << x0[t][j] << "\t";
			out << endl;
		}
		out << endl;

		double r = sqrt(reciever_x * reciever_x + reciever_y * reciever_y);
		double B2D = Get_solution(t, r, reciever_z, -reciever_y, true);

	fout << setprecision(16) << time[t] << "\t" << B2D << endl;		//ищем B общего поля в приемнике
	}
	//------------------------------------------------------------------------------------------

	//-------------------Инициализация решения на 3-м и посл. слоях-----------------------------
	//for (int t = 2; t < time_number; t++)
	//{
	//	Clear_and_resize(); //очищаем вектора для каждого временного слоя

	//	//считаем коэффициенты по времени
	//	t0 = time[t] - time[t - 1]; //t0
	//	t1 = time[t] - time[t - 2]; //t
	//	t2 = time[t - 1] - time[t - 2]; //t1
	//	nu0 = (t1 + t0) / t1 / t0;
	//	nu1 = t1 / t2 / t0;
	//	nu2 = t0 / t2 / t1;

	//	for (int i = 0; i < Kel; i++)
	//	{
	//		Locals(i, time[t]);
	//		for (int i = 0; i < rasm; i++)
	//			for (int j = 0; j < rasm; j++)
	//				A_loc[i][j] = G_loc[i][j] + M_loc[i][j] * nu0 + M_loc_g[i][j];

	//		Assemble_Locals(i); //сборка левой части

	//		//сборка вектора правой части F = b_loc(j)-nu2*M*x2+nu1*M*x1
	//		for (int j = 0; j < rasm; j++)
	//			temp[j] = x0[t-2][elems[i].num[j]];
	//		multiply(M_loc, temp, M2);//M*x2
	//		for (int j = 0; j < rasm; j++)
	//			temp[j] = x0[t-1][elems[i].num[j]];
	//		multiply(M_loc, temp, M1);//M*x1
	//		for (int j = 0; j < rasm; j++) //зачем этот цикл???
	//		{
	//			M1[j] = M1[j] * nu1;
	//			M2[j] = M2[j] * nu2;
	//		}
	//		for (int j = 0; j < rasm; j++)
	//			b[elems[i].num[j]] += b_loc[j] - M2[j] + M1[j];
	//	}

	//	Get_KR1(time[t]);

	//	LOS_LU(t);
	//	f << "t = " << time[t] << endl;
	//	for (int i = 0; i < N; i++)
	//		f << nodes[i].r << "\t" << nodes[i].z << "\t" << x0[t][i] << endl;
	//	f << endl;

	//	out << "t = " << time[t] << endl;
	//	for (int i = 1; i <= N / num_r; i++)
	//	{
	//		for (int j = N - num_r * i; j < N - num_r * (i - 1); j++)
	//			out << x0[t][j] << "\t";
	//		out << endl;
	//	}
	//	out << endl;

	////	double B2D = Get_solution(t, r, reciever_z, -reciever_y, true);
	////	fout << setprecision(16) << time[t] << "\t" << B2D << endl;
	//}
	out.close();
	f.close();
//	fout.close();
}

void KURS::Get_KR1(double t) // учет первых краевых
{
		for (int j = 0; j < KR1.size(); j++)
		{
			int node_id = KR1[j];
			di[node_id] = 1;
			b[node_id] = 0.;
			for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
				al[k] = 0;
			for (int k = 0; k < ja.size(); k++)
				if (ja[k] == node_id)
					au[k] = 0;
		}
}
//----------------------------------------------------------------------

////-------------блок функций и процедур для 3D задачи--------------------

double KURS::Get_solution(int time_layer, double r, double z, double variable, bool isB)
{
		bool finded = false;
		int fe_sol = 0;
		for (int i = 0; i < num_r*num_z && !finded; i++) //сделать разные переменные для N!!!
		{
			if (r >= nodes[elems[i].num[0]].r && r <= nodes[elems[i].num[1]].r &&
				z >= nodes[elems[i].num[0]].z && z <= nodes[elems[i].num[2]].z)
			{
				finded = true;
				fe_sol = i;
			}
		}

		// Если нашли, то решение будет линейной комбинацией базисных функций на соответствующие веса
		// Вычисление шага
		double hr = fabs(nodes[elems[fe_sol].num[1]].r - nodes[elems[fe_sol].num[0]].r);
		//cout << hr << " ";
		double hz = fabs(nodes[elems[fe_sol].num[2]].z - nodes[elems[fe_sol].num[0]].z);
		//cout << hz << " ";
		double A0 = 0.0, A1 = 0.0, result = 0.0;

		if (isB == false)
		{
			// Находим линейные одномерные функции
			double X1 = (nodes[elems[fe_sol].num[1]].r - r) / hr;
			double X2 = (r - nodes[elems[fe_sol].num[0]].r) / hr;
			double Y1 = (nodes[elems[fe_sol].num[2]].z - z) / hz;
			double Y2 = (z - nodes[elems[fe_sol].num[0]].z) / hz;

			// Находим значение билинейных базисных функций
			double psi[4];
			psi[0] = X1 * Y1;
			psi[1] = X2 * Y1;
			psi[2] = X1 * Y2;
			psi[3] = X2 * Y2;

			// Линейная комбинация базисных функций на веса
			for (int i = 0; i < 4; i++)
			{
				A0 += x0[time_layer - 1][elems[fe_sol].num[i]] * psi[i]; //на предыдущем слое
				A1 += x0[time_layer][elems[fe_sol].num[i]] * psi[i]; //на текущем слое
			}
			A0 = variable / r * A0;
			A1 = variable / r * A1;
			result = (A0 - A1) / t0; //t0 = t[i]-t[i-1], то есть шаг по времени
		}
		else
		{
			// Находим линейные одномерные функции
			double X1 = (nodes[elems[fe_sol].num[1]].r - r) / hr;
			double X2 = (r - nodes[elems[fe_sol].num[0]].r) / hr;
			double Y1 = (nodes[elems[fe_sol].num[2]].z - z) / hz;
			double Y2 = (z - nodes[elems[fe_sol].num[0]].z) / hz;

			// Находим значение билинейных базисных функций
			double psi_der[8];
			psi_der[0] = X1 * Y1;
			//cout << psi_der[0] << " ";
			psi_der[1] = X2 * Y1;
		//	cout << psi_der[1] << " ";
			psi_der[2] = X1 * Y2;
			//cout << psi_der[2] << " ";
			psi_der[3] = X2 * Y2;
		//	cout << psi_der[3] << " ";
			psi_der[4] = -Y1 / hr;
			//cout << psi_der[4] << " ";
			psi_der[5] = Y1 / hr;
			//cout << psi_der[5] << " ";
			psi_der[6] = -Y2 / hr;
			//cout << psi_der[6] << " ";
			psi_der[7] = Y2 / hr;
			//cout << psi_der[7] << " ";

			// Линейная комбинация базисных функций на веса
			for (int i = 0; i < 4; i++)
			{
				A0 = A0 + x0[time_layer][elems[fe_sol].num[i]] * psi_der[i]; //считаем А
				A1 += x0[time_layer][elems[fe_sol].num[i]] * psi_der[i+4]; //Производная по А
				//cout << x0[time_layer][elems[fe_sol].num[i]] << " " << A1 << " ";
			}
			/*A1 = x0[time_layer][elems[fe_sol].num[0]] * psi_der[4] * (nodes[elems[fe_sol].num[1]].r / r - 2) +
				x0[time_layer][elems[fe_sol].num[1]] * psi_der[4] * (2- nodes[elems[fe_sol].num[0]].r/r) +
				x0[time_layer][elems[fe_sol].num[2]] * psi_der[6] * (nodes[elems[fe_sol].num[1]].r / r - 2) +
				x0[time_layer][elems[fe_sol].num[3]] * psi_der[6] * (2 - nodes[elems[fe_sol].num[0]].r / r);*/
			result = A0/r + A1; 
		}
	//	cout << result << endl;
		return result;
}

void KURS::Middles()
{
	ofstream fout("xyz_coords.txt");
	double middle = 0.0, diff = 0.0;
	fout << edges.size() << endl;
	for (int i = 0; i < edges.size(); i++)
	{
		diff = edges[i].end.x - edges[i].beg.x;
		middle = edges[i].beg.x + diff / 2.;
		if (diff != 0)
		{
			fout << middle << " " << edges[i].beg.y << " " << edges[i].beg.z << endl;
			edges[i].middle.x = middle;
			edges[i].middle.y = edges[i].beg.y;
			edges[i].middle.z = edges[i].beg.z;
		}
		else
		{
			diff = edges[i].end.y - edges[i].beg.y;
			middle = edges[i].beg.y + diff / 2.;
			if (diff != 0)
			{
				fout << edges[i].beg.x << " " << middle << " " << edges[i].beg.z << endl;
				edges[i].middle.x = edges[i].beg.x;
				edges[i].middle.y = middle;
				edges[i].middle.z = edges[i].beg.z;
			}
			else
			{
				diff = edges[i].end.z - edges[i].beg.z;
				middle = edges[i].beg.z + diff / 2.;
				fout << edges[i].beg.x << " " << edges[i].beg.y << " " << middle << endl;
				edges[i].middle.x = edges[i].beg.x;
				edges[i].middle.y = edges[i].beg.y;
				edges[i].middle.z = middle;
			}
		}
	}
	fout.close();
}

void KURS::Generate_Portrait_3D() // генерация портрета
{
	N = num_edges;
	Kel = num_elems;
	ia.resize(N + 1);
	ja.resize(rasm * rasm * Kel);
	vector<int> t1(rasm * rasm * Kel);
	vector<int> t2(rasm * rasm * Kel);
	vector<int> beg(N);
	int s = 0;
	for (int i = 0; i < N; i++)
		beg[i] = 0;
	for (int el = 0; el < Kel; el++)
	{
		for (int i = 0; i < rasm; i++)
		{
			int k = obj_elems[el].num[i];
			for (int j = i + 1; j < rasm; j++)
			{
				int ind1 = k;
				int ind2 = obj_elems[el].num[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				int iaddr = beg[ind2];
				if (iaddr == 0)
				{
					s++;
					beg[ind2] = s;
					t1[s] = ind1;
					t2[s] = 0;
				}
				else
				{
					while (t1[iaddr] < ind1 && t2[iaddr] > 0)
						iaddr = t2[iaddr];
					if (t1[iaddr] > ind1)
					{
						s++;
						t1[s] = t1[iaddr];
						t2[s] = t2[iaddr];
						t1[iaddr] = ind1;
						t2[iaddr] = s;
					}
					else if (t1[iaddr] < ind1)
					{
						s++;
						t2[iaddr] = s;
						t1[s] = ind1;
						t2[s] = 0;
					}
				}
			}
		}
	}

	ia[0] = 0;
	for (int i = 0; i < N; i++)
	{
		ia[i + 1] = ia[i];
		int a = beg[i];
		while (a != 0)
		{
			ja[ia[i + 1]] = t1[a];
			ia[i + 1]++;
			a = t2[a];
		}
	}
}

void KURS::Make_grid(bool go_to_z_up)
{
	//строим сетку по x
	num_x = 0; //обнуляем число элементов по x 
	double h = h_obj_x / d_x;
	bool flag = true;
	double usel = x_obj_left;

	ax.push_back(x_obj_left);
	for (int i = 0; usel > x_left; i++) //идем влево от источника
	{
		flag = true;
		h = h * d_x;
		usel = ax[i] - h;
		if (usel < x_left)
		{
			flag = false;
			usel = x_left;
			if (abs(ax[i] - usel) < abs(ax[i] - ax[i - 1]))
			{
				ax[i] = usel;
				h = abs(ax[i] - ax[i - 1]) / d_x;
				i--;
			}
			else
				ax.push_back(usel);
		}
		if (flag == true)
			ax.push_back(usel);
	}

	reverse(ax.begin(), ax.end());

	usel = x_obj_left;
	h = h_obj_x / d_obj_x;
	for (int i = ax.size()-1; usel < x_obj_right; i++) //идем вправо от начала координат
	{
		flag = true;
		h = h * d_obj_x;
		usel = ax[i] + h;
		if (usel > x_obj_right)
		{
			usel = x_obj_right;
			flag = false;
			if (abs(ax[i] - usel) < abs(ax[i] - ax[i - 1]))
			{
				ax[i] = usel;
				h = abs(ax[i] - ax[i - 1]) / d_obj_x;
				i--;
			}
			else
				ax.push_back(usel);
		}
		if (flag == true)
			ax.push_back(usel);
	}

	h = h / d_x;
	for (int i = ax.size()-1; usel < x_right; i++) //идем вправо от объекта
	{
		flag = true;
		h = h * d_x;
		usel = ax[i] + h;
		if (usel > x_right)
		{
			usel = x_right;
			flag = false;
			if (abs(ax[i] - usel) < abs(ax[i] - ax[i - 1]))
			{
				ax[i] = usel;
				h = abs(ax[i] - ax[i - 1]) / d_x;
				i--;
			}
			else
				ax.push_back(usel);
		}
		if (flag == true)
			ax.push_back(usel);
	}
	num_x = ax.size();

	for (int j = 0; j < num_x; j++)
		x_s.insert(ax[j]);

	//строим сетку по y 
	num_y = 0; //обнуляем число элементов по y
	ay.push_back(y_obj_left);
	usel = y_obj_left;
	h = h_obj_y / d_y;

	for (int i = 0; usel > y_left; i++) //идем влево от источника
	{
		flag = true;
		h = h * d_y;
		usel = ay[i] - h;
		if (usel < y_left)
		{
			flag = false;
			usel = y_left;
			if (abs(ay[i] - usel) < abs(ay[i] - ay[i - 1]))
			{
				ay[i] = usel;
				h = abs(ay[i] - ay[i - 1]) / d_y;
				i--;
			}
			else
				ay.push_back(usel);
		}
		if (flag == true)
			ay.push_back(usel);
	}

	reverse(ay.begin(), ay.end());

	h = h_obj_y / d_obj_y;
	usel = y_obj_left;

	for (int i = ay.size()-1; usel < y_obj_right; i++) //идем от начала координат
	{
		flag = true;
		h = h * d_obj_y;
		usel = ay[i] + h;
		if (usel > y_obj_right)
		{
			usel = y_obj_right;
			flag = false;
			if (abs(ay[i] - usel) < abs(ay[i] - ay[i - 1]))
			{
				ay[i] = usel;
				h = abs(ay[i] - ay[i - 1]) / d_obj_y;
				i--;
			}
			else
				ay.push_back(usel);
		}
		if (flag == true)
			ay.push_back(usel);
	}

	h = h / d_y;
	for (int i = ay.size()-1; usel < y_right; i++) //идем от объекта
	{
		flag = true;
		h = h * d_y;
		usel = ay[i] + h;
		if (usel > y_right)
		{
			usel = y_right;
			flag = false;
			if (abs(ay[i] - usel) < abs(ay[i] - ay[i - 1]))
			{
				ay[i] = usel;
				h = abs(ay[i] - ay[i - 1]) / d_y;
				i--;
			}
			else
				ay.push_back(usel);
		}
		if (flag == true)
			ay.push_back(usel);
	}

	num_y = ay.size();

	//обработка по y
	for (int j = 0; j < num_y; j++)
		y_s.insert(ay[j]);

	if (go_to_z_up == true) //для уплотнения сетки по z дальше от источника
	{
		//строим сетку по z 
		num_obj_z = 0; //обнуляем число элементов по z
		a_obj_z.push_back(z_obj_bottom);
		usel = z_obj_bottom;
		h = h_obj_z / d_obj_z;
		for (int i = 0; usel < z_obj_top; i++) //идем вверх от начала координат
		{
			flag = true;
			h = h * d_obj_z;
			usel = a_obj_z[i] + h;

			if (usel > z_obj_top)
			{
				usel = z_obj_top;
				flag = false;
				if (abs(a_obj_z[i] - usel) < abs(a_obj_z[i] - a_obj_z[i - 1]))
				{
					a_obj_z[i] = usel;
					h = abs(a_obj_z[i] - a_obj_z[i - 1]) / d_obj_z;
					i--;
				}
				else
					a_obj_z.push_back(usel);
			}
			if (flag == true)
				a_obj_z.push_back(usel);
		}
	}
	else
	{
		//для уплотнения сетки по ближе к источнику
		num_obj_z = 0; //обнуляем число элементов по z
		a_obj_z.push_back(z_obj_top);
		usel = z_obj_top;
		h = h_obj_z / d_obj_z;
		for (int i = 0; usel > z_obj_bottom; i++) //идем вверх от начала координат
		{
			flag = true;
			h = h * d_obj_z;
			usel = a_obj_z[i] - h;
			if (usel < z_obj_bottom)
			{
				usel = z_obj_bottom;
				flag = false;
				if (abs(a_obj_z[i] - usel) < abs(a_obj_z[i] - a_obj_z[i - 1]))
				{
					a_obj_z[i] = usel;
					h = abs(a_obj_z[i] - a_obj_z[i - 1]) / d_obj_z;
					i--;
				}
				else
					a_obj_z.push_back(usel);
			}
			if (flag == true)
				a_obj_z.push_back(usel);
		}

		h = h / d_sr_z;
		for (int i = a_obj_z.size()-1; usel > z_sr_left; i++) //идем вверх от начала координат
		{
			flag = true;
			h = h * d_sr_z;
			usel = a_obj_z[i] - h;
			for (int j = 0; j < n_layers; j++) //проеряем значения в окрестностях границ слоев
			if (layers[j].y0 > usel && a_obj_z[i] > layers[j].y0)
			{
				usel = layers[j].y0;
				flag = false;
				if (abs(a_obj_z[i] - usel) < abs(a_obj_z[i] - a_obj_z[i - 1]))
				{
					a_obj_z[i] = usel;
					h = abs(a_obj_z[i] - a_obj_z[i - 1]) / d_sr_z;
					i--;
				}
				else
					a_obj_z.push_back(usel);
			}
			if (flag == true)
				a_obj_z.push_back(usel);
		}
		reverse(a_obj_z.begin(), a_obj_z.end());

		h = h_obj_z / d_sr_z;
		for (int i = a_obj_z.size() - 1; usel < z_sr_right; i++) //идем от начала координат
		{
			flag = true;
			h = h * d_sr_z;
			usel = a_obj_z[i] + h;
			for (int j = 0; j < n_layers; j++)
			if (layers[j].y1 < usel && a_obj_z[i] < layers[j].y1)
			{
				usel = layers[j].y1;
				flag = false;
				if (abs(a_obj_z[i] - usel) < abs(a_obj_z[i] - a_obj_z[i - 1]))
				{
					a_obj_z[i] = usel;
					h = abs(a_obj_z[i] - a_obj_z[i - 1]) / d_sr_z;
					i--;
				}
				else
					a_obj_z.push_back(usel);
			}
			if (flag == true)
				a_obj_z.push_back(usel);
		}
	}


	num_obj_z = a_obj_z.size();

	//обработка по z
	for (int j = 0; j < num_obj_z; j++)
		z_obj_s.insert(a_obj_z[j]);

}

void KURS::Make_edges() //Формирование ребер
{
	num_edges = num_obj_z * (num_x - 1) * num_y + num_obj_z * (num_y - 1) * num_x + (num_obj_z - 1) * num_y * num_x;
	edges.resize(num_edges);
	int l = 0;

	// ребра параллельные х
	for (int i = 0; i < num_obj_z; i++)
		for (int j = 0; j < num_x - 1; j++)
			for (int k = 0; k < num_y; k++)
			{
				edges[l].beg.x = ax[j];
				edges[l].beg.y = ay[k];
				edges[l].beg.z = a_obj_z[i];
				edges[l].end.x = ax[j + 1];
				edges[l].end.y = ay[k];
				edges[l].end.z = a_obj_z[i];
				l++;
			}
	// ребра параллельные y
	for (int i = 0; i < num_obj_z; i++)
		for (int j = 0; j < num_y - 1; j++)
			for (int k = 0; k < num_x; k++)
			{
				edges[l].beg.x = ax[k];
				edges[l].beg.y = ay[j];
				edges[l].beg.z = a_obj_z[i];
				edges[l].end.x = ax[k];
				edges[l].end.y = ay[j + 1];
				edges[l].end.z = a_obj_z[i];
				l++;
			}
	// ребра параллельные z
	for (int i = 0; i < num_obj_z - 1; i++)
		for (int k = 0; k < num_y; k++)
			for (int j = 0; j < num_x; j++)
			{
				edges[l].beg.x = ax[j];
				edges[l].beg.y = ay[k];
				edges[l].beg.z = a_obj_z[i];
				edges[l].end.x = ax[j];
				edges[l].end.y = ay[k];
				edges[l].end.z = a_obj_z[i + 1];
				l++;
			}

}

int KURS::Get_Num_Layer_3D(double z0, double z1)
{
	for (int i = 0; i < n_layers; i++)
		if (layers[i].y0 <= z0 && layers[i].y1 >= z1)
			return i;
	return 0;
}

void KURS::Make_obj_elems()
{
	int shift_y = num_obj_z * num_y * (num_x - 1),
		shift_z = shift_y + num_obj_z * num_x * (num_y - 1), l = 0;

	num_elems = (num_obj_z - 1) * (num_x - 1) * (num_y - 1);
	obj_elems.resize(num_elems);

	for (int i = 0; i < num_obj_z - 1; i++)
		for (int j = 0; j < num_x - 1; j++)
			for (int k = 0; k < num_y - 1; k++)
			{
				obj_elems[l].num.resize(12);
				obj_elems[l].num[0] = i * num_y * (num_x - 1) + j * num_y + k;
				obj_elems[l].num[1] = i * num_y * (num_x - 1) + j * num_y + k + 1;
				obj_elems[l].num[2] = (i + 1) * num_y * (num_x - 1) + j * num_y + k;
				obj_elems[l].num[3] = (i + 1) * num_y * (num_x - 1) + j * num_y + k + 1;
				obj_elems[l].num[4] = shift_y + i * num_x * (num_y - 1) + j + k * num_x;
				obj_elems[l].num[5] = shift_y + i * num_x * (num_y - 1) + j + 1 + k * num_x;
				obj_elems[l].num[6] = shift_y + (i + 1) * num_x * (num_y - 1) + j + k * num_x;
				obj_elems[l].num[7] = shift_y + (i + 1) * num_x * (num_y - 1) + j + 1 + k * num_x;
				obj_elems[l].num[8] = shift_z + i * num_y * num_x + k * num_x + j;
				obj_elems[l].num[9] = shift_z + i * num_y * num_x + k * num_x + 1 + j;
				obj_elems[l].num[10] = shift_z + i * num_y * num_x + (k + 1) * num_x + j;
				obj_elems[l].num[11] = shift_z + i * num_y * num_x + (k + 1) * num_x + 1 + j;

				if(edges[obj_elems[l].num[0]].end.x <= x_obj_right &&
					edges[obj_elems[l].num[0]].beg.x >= x_obj_left &&
					edges[obj_elems[l].num[4]].end.y <= y_obj_right &&
					edges[obj_elems[l].num[4]].beg.y >= y_obj_left &&
					edges[obj_elems[l].num[8]].end.z <= z_obj_top &&
					edges[obj_elems[l].num[8]].beg.z >= z_obj_bottom)
					obj_elems[l].mater = Nmat; //следующий за слоями 2d сетки
				else
				{
					unsigned int num_area = Get_Num_Layer_3D(edges[obj_elems[l].num[8]].beg.z, edges[obj_elems[l].num[8]].end.z);
					obj_elems[l].mater = layers[num_area].n_mat;
				}
					/*obj_elems[l].mater = 1;
				else if (edges[obj_elems[l].num[8]].end.z <= z_obj_bottom)
					obj_elems[l].mater = 2;
				else 
					obj_elems[l].mater = 0;*/
				l++;
			}
}

void KURS::Make_centers()
{
	//формируем точки - середины элементов, в которых будем искать значения
	centers.resize(num_elems);
	for (int i = 0; i < num_elems; i++)
	{
		element el = obj_elems[i];
		hx = edges[el.num[0]].end.x - edges[el.num[0]].beg.x;
		hy = edges[el.num[4]].end.y - edges[el.num[4]].beg.y;
		hz_o = edges[el.num[8]].end.z - edges[el.num[8]].beg.z;
		centers[i].x = edges[el.num[0]].beg.x + hx / 2.;
		centers[i].y = edges[el.num[4]].beg.y + hy / 2.;
		centers[i].z = edges[el.num[8]].beg.z + hz_o / 2.;
	}
}

void KURS::Generate_3D_Net(string net_file, bool go_to_z_up)
{
	ifstream fin(net_file + ".txt");
	fin >> x_obj_left >> x_obj_right >> h_obj_x >> d_obj_x >>
		y_obj_left >> y_obj_right >> h_obj_y >> d_obj_y >>
		z_obj_bottom >> z_obj_top >> h_obj_z >> d_obj_z >>
		lay_obj >> 
		x_left >> x_right >> d_x >>
		y_left >> y_right >> d_y >>
		z_sr_left >> z_sr_right >> d_sr_z; //считываем границы области, шаги по x,y,z, коэффициент растяжения, слой снизу, в котором лежит объект

	Make_grid(go_to_z_up);
	Make_edges();
	Make_obj_elems();
	Make_centers();

	// параллельно х
	// верх - низ
	for (int k = 0; k < num_y; k++)
	{
		for (int j = 0; j < num_x - 1; j++)
		{
			KR1_3D.push_back((num_obj_z - 1) * num_y * (num_x - 1) + j * num_y + k);
			KR1_3D.push_back(j * num_y + k);
		}
	}

	// лево - право
	for (int i = 1; i < num_obj_z - 1; i++)
	{
		for (int j = 0; j < num_x - 1; j++)
		{
			KR1_3D.push_back(i * num_y * (num_x - 1) + j * num_y);
			KR1_3D.push_back(i * num_y * (num_x - 1) + (j + 1) * num_y - 1);
		}
	}

	// параллельно у
	// верх - низ
	int shift_y = num_obj_z * num_y * (num_x - 1),
		shift_z = shift_y + num_obj_z * num_x * (num_y - 1), l = 0;

	for (int k = 0; k < num_y - 1; k++)
	{
		for (int j = 0; j < num_x; j++)
		{
			KR1_3D.push_back(shift_y + (num_obj_z - 1) * num_x * (num_y - 1) + k * num_x + j);
			KR1_3D.push_back(shift_y + k * num_x + j);
		}
	}
	// перед - зад
	for (int i = 1; i < num_obj_z - 1; i++)
	{
		for (int k = 0; k < num_y - 1; k++)
		{
			KR1_3D.push_back(shift_y + i * num_x * (num_y - 1) + k * num_x);
			KR1_3D.push_back(shift_y + i * num_x * (num_y - 1) + (k + 1) * num_x - 1);
		}
	}

	// параллельно z
	// лево - право
	for (int i = 0; i < num_obj_z - 1; i++)
	{
		for (int j = 0; j < num_x; j++)
		{
			KR1_3D.push_back(shift_z + i * num_y * num_x + j);
			KR1_3D.push_back(shift_z + i * num_y * num_x + num_x * (num_y - 1) + j);
		}
	}
	// перед - зад
	for (int i = 0; i < num_obj_z - 1; i++)
	{
		for (int k = 1; k < num_y - 1; k++)
		{
			KR1_3D.push_back(shift_z + i * num_y * num_x + k * num_x);
			KR1_3D.push_back(shift_z + i * num_y * num_x + (k + 1) * num_x - 1);
		}
	}
}

double KURS::Get_solution_3D(double t, double x, double y, double z, int i, bool B)
{

		element el = obj_elems[i];
		hx = edges[el.num[0]].end.x - edges[el.num[0]].beg.x;
		hy = edges[el.num[4]].end.y - edges[el.num[4]].beg.y;
		hz_o = edges[el.num[8]].end.z - edges[el.num[8]].beg.z;

		double result_x = 0.0, result_y = 0.0, result_z = 0.0, modul;

		// Находим линейные одномерные функции
		double X1 = (edges[el.num[0]].end.x - x) / hx;
		double X2 = (x - edges[el.num[0]].beg.x) / hx;
		double Y1 = (edges[el.num[1]].beg.y - y) / hy;
		double Y2 = (y - edges[el.num[0]].beg.y) / hy;
		double Z1 = (edges[el.num[3]].beg.z - z) / hz_o;
		double Z2 = (z - edges[el.num[0]].beg.z) / hz_o;

		if (B == false)
		{
			// Находим значение триилинейных базисных функций
			double psi[12];
			psi[0] = Y1 * Z1; 
			psi[1] = Y2 * Z1;
			psi[2] = Y1 * Z2;
			psi[3] = Y2 * Z2;

			psi[4] = X1 * Z1;
			psi[5] = X2 * Z1;
			psi[6] = X1 * Z2;
			psi[7] = X2 * Z2;

			psi[8] = X1 * Y1;
			psi[9] = X2 * Y1;
			psi[10] = X1 * Y2;
			psi[11] = X2 * Y2;

			// Линейная комбинация базисных функций на веса
			for (int i = 0; i < 4; i++)
				result_x += x0_3D[t][el.num[i]] * psi[i];

			for (int i = 4; i < 8; i++)
				result_y += x0_3D[t][el.num[i]] * psi[i];

			for (int i = 8; i < 12; i++)
				result_z += x0_3D[t][el.num[i]] * psi[i];

			modul = sqrt(result_x * result_x + result_y * result_y + result_z * result_z);
		}
		else
		{
			double psi[8]; 
			psi[0] = -Z1/hy; //dx/dy
			psi[1] =  Z1/hy;
			psi[2] = -Z2/hy;
			psi[3] =  Z2/hy;

			psi[4] = -Z1/hx; //dy/dx
			psi[5] = Z1/hx;
			psi[6] = -Z2/hx;
			psi[7] = Z2/hx;

			for (int i = 0; i < 4; i++)
				result_x += x0_3D[t][el.num[i]] * psi[i];

			for (int i = 4; i < 8; i++)
				result_y += x0_3D[t][el.num[i]] * psi[i];

			modul = result_y - result_x;
		}
			return modul;
}

void KURS::Locals_3D(int el_id, int time_layer) // получение локальной матрицы А
{
	element el = obj_elems[el_id];
	hx = edges[el.num[0]].end.x - edges[el.num[0]].beg.x;
	hy = edges[el.num[4]].end.y - edges[el.num[4]].beg.y;
	hz_o = edges[el.num[8]].end.z - edges[el.num[8]].beg.z;
	lambda = materials[el.mater].lambda;
	sigma = materials[el.mater].sigma;
	Get_G_3D();
	Get_M_3D();
	Get_b_3D(time_layer, el_id);
}

void KURS::Get_G_3D() // получение локальной G
{
	double a1 = (lambda * hx * hy) / (6 * hz_o),
		a2 = (lambda * hx * hz_o) / (6 * hy),
		a3 = (lambda * hz_o * hy) / (6 * hx),
		a4 = -(lambda * hz_o) / 6.,
		a5 = (lambda * hy) / 6., 
		a6 = -(lambda * hx) / 6.;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		     G_loc[i][j] = a1 * G1[i][j] + a2 * G2[i][j];
	
	for (int i = 4; i < 8; i++)
		for (int j = 0; j < 4; j++)
			G_loc[j][i] = G_loc[i][j] = a4 * G2[i-4][j];

	for (int i = 8; i < 12; i++)
		for (int j = 0; j < 4; j++)
			G_loc[i][j] = a5 * G3T[i - 8][j];

	for (int i = 4; i < 8; i++)
		for (int j = 4; j < 8; j++)
			G_loc[i][j] = a1 * G1[i - 4][j-4] + a3*G2[i-4][j-4];

	for (int i = 8; i < 12; i++)
		for (int j = 4; j < 8; j++)
			G_loc[j][i] = G_loc[i][j] = a6 * G1[i - 8][j-4];

	for (int i = 0; i < 4; i++)
		for (int j = 8; j < 12; j++)
			G_loc[i][j] = a5 * G3[i][j-8];

	for (int i = 8; i < 12; i++)
		for (int j = 8; j < 12; j++)
			G_loc[i][j] = a2 * G1[i - 8][j - 8] + a3 * G2[i - 8][j - 8];
}

void KURS::Get_M_3D() // получение локальной М
{
	double koef = hx * hy * hz_o / 36.;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			M_loc[i][j] = M_loc[j][i] = koef * M[i][j];
	
	for (int i = 4; i < 8; i++)
		for (int j = 4; j < 8; j++)
			M_loc[i][j] = M_loc[j][i] = koef * M[i-4][j-4];

	for (int i = 8; i < 12; i++)
		for (int j = 8; j < 12; j++)
			M_loc[i][j] = M_loc[j][i] = koef * M[i-8][j-8];
}

void KURS::Get_b_3D(int time_layer, int el_id) // получение локального b
{
	F.clear();
	F.resize(12);
	element el = obj_elems[el_id];
	if (el.mater == Nmat && time_layer != 0) //для области объекта
	{
		double x, y, r, z;
		double koef_sigma = sigma - materials[layers[lay_obj].n_mat].sigma;
		for (int i = 0; i < 8; i++)
		{
			x = edges[el.num[i]].middle.x;
			y = edges[el.num[i]].middle.y;
			r = sqrt(x * x + y * y);
			z = edges[el.num[i]].middle.z;
			if (i < 4)
				F[i] = koef_sigma * Get_solution(time_layer, r, z, -y, false);
			else
				F[i] = koef_sigma * Get_solution(time_layer, r, z, x, false);
		}
		for (int i = 8; i < 11; i++)
			F[i] = 0.;
	}
	//b = C*F, где C = М без сигмы
	multiply(M_loc, F, b_loc);
}

void KURS::Assemble_Locals_3D(int el_id) // внесение локальных A, b  в глобальную СЛАУ
{
	vector<int> L = obj_elems[el_id].num;
	int k = obj_elems[el_id].num.size(); // размерность локальной матрицы
	for (int i = 0; i < k; i++)
		di[L[i]] += A_loc[i][i];

	for (int i = 0; i < rasm; i++)
		for (int j = 0; j < i; j++)
			for (int k = ia[L[i]]; k < ia[L[i] + 1]; k++)
				if (ja[k] == L[j])
				{
					al[k] += A_loc[i][j];
					au[k] += A_loc[j][i];
					k++;
					break;
				}
}

void KURS::Find_in_object()
{
//	cout << endl << num_x << " " << num_y << " " << num_obj_z << endl;
	N = num_edges;

	x0_3D.resize(time_number);
	for (int i = 0; i < time_number; i++)
		x0_3D[i].resize(N);
	temp.resize(N);

	int center_y = 41;
	/*int center_y = 82;*/
	//-------Инициализация решения на 0-м слое как решение стационарной задачи----------------
	Clear_and_resize();//очищаем вектора для каждого временного слоя
	ofstream ff("chekelems.txt");

	for (int i = 0; i < num_elems; i++) 
	{
		ff << sigma << " ";
		ff << edges[obj_elems[i].num[8]].beg.z << " " << edges[obj_elems[i].num[8]].end.z << endl;
		Locals_3D(i, 0);
		/*if (obj_elems[i].mater == Nmat)
			sigma = materials[layers[lay_obj].n_mat].sigma;*/
		for (int i = 0; i < rasm; i++)
			for (int j = 0; j < rasm; j++)
				A_loc[i][j] = G_loc[i][j] + sigma*M_loc[i][j];

		Assemble_Locals_3D(i); //сборка левой части

		for (int j = 0; j < rasm; j++) 		//сборка вектора правой части
			b[obj_elems[i].num[j]] += b_loc[j];
	}
	Get_KR1_3D();
	//diag_preconditioning(x0_3D[0]);
	LOS_LU_3D(0);
//	z.clear();
//	L.clear();
//	D.clear();
//	U.clear();
//	z.resize(N);
//	L.resize(ia[N]);
//	D.resize(N);
//	U.resize(ia[N]);
//	//LOS_LU_3D(t);
////	diag_preconditioning(x0_3D[t]);
//	FactLU(L, U, D);
//	Direct(L, D, z, b);
//	Reverse(U, x0_3D[0], z);
	ofstream out("result_object.txt");
	//for (int i = 0; i < N; i++) //выводим результат в ребрах
	//	out << x0_3D[0][i] << endl;
	//out << endl;

	ofstream f("y_srez.txt");
	bool finded = false;
	int reciever_elem = 0;
	double otv;

	f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
	ofstream xout("result_xz.txt");

	out << "t = " << time[0] << endl;
	for (int i = 0; i < num_elems && !finded; i++)
	{
		if (20 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			20 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			-20 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		{
			finded = true;
			reciever_elem = i;
		}
	}

	otv = Get_solution_3D(0,20, 20, -20, reciever_elem, false);
	out << "20 20 -20" << "\t" << otv << endl;
	finded = false;
	for (int i = 0; i < num_elems && !finded; i++)
	{
		if (0 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			0 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			0 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		{
			finded = true;
			reciever_elem = i;
		}
	}

	//f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
	//ofstream xout("result_xz.txt");
	otv = Get_solution_3D(0, 0, 0, 0, reciever_elem, false);
	out << "0 0 0" << "\t" << otv << endl;
	finded = false;
	for (int i = 0; i < num_elems && !finded; i++)
	{
		if (50 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			50 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			-50 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		{
			finded = true;
			reciever_elem = i;
		}
	}

	//f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
	//ofstream xout("result_xz.txt");
	otv = Get_solution_3D(0, 50, 50, -50, reciever_elem, false);
	out << "50 50 -50" << "\t" << otv << endl;

	//		otvet.push_back(otv);
	//for (int i = 0; i < num_elems; i++)
	//	if (centers[i].y == centers[center_y].y) //срез по y
	//	{
	//		otv = Get_solution_3D(0, centers[i].x, centers[i].y, centers[i].z, i, false);
	//		otvet.push_back(otv);
	//		if (centers[i].x > -30 && centers[i].x < 60 && centers[i].z < 20 && centers[i].z > -60)
	//		out << centers[i].x << "\t" << centers[i].z << "\t"  << otv << endl;
	//	}

	//for (int i = 0; i < otvet.size(); i++) //выводим результат в ребрах
	//   out << centers[i].x << "\t" << centers[i].z << "\t" << otvet[i] << endl;
	//out << endl;

	for (int i = 1; i <= otvet.size() / (num_x - 1); i++) // вывод в файл среза по y
	{
		for (int j = otvet.size() - (num_x - 1) * i; j < otvet.size() - (num_x - 1) * (i - 1); j++)
			f << otvet[j] << "\t";
		f << endl;
	}
	f << endl;

	otvet.clear();
	//------------------------------------------------------------------------------------------
	ofstream fout("in_reciever.txt");
	double r = sqrt(reciever_x * reciever_x + reciever_y * reciever_y);
	//-----------------------Инициализация решения на 1-м слое----------------------------------
	for (int t = 1; t < 2; t++)
	{
	//int t = 1;
	   Clear_and_resize();//очищаем вектора для каждого временного слоя

	   //считаем коэффициенты по времени для двухслойной неявной схемы
	   t0 = time[t] - time[t-1]; //t0
	   nu0 = 1 / t0;
	 // / cout << nu0 << endl;

	   for (int i = 0; i < num_elems; i++)
	   {
		  Locals_3D(i, t);
		  /*	if (obj_elems[i].mater == Nmat)
				  sigma = materials[layers[lay_obj].n_mat].sigma;*/
		  for (int i = 0; i < rasm; i++)
			 for (int j = 0; j < rasm; j++)
				A_loc[i][j] = G_loc[i][j] + sigma * M_loc[i][j] * nu0;

		  Assemble_Locals_3D(i); //сборка левой части

		  //сборка вектора правой части
		  for (int j = 0; j < rasm; j++)
			 temp[j] = x0_3D[t-1][obj_elems[i].num[j]];
		  multiply(M_loc, temp, M2);//M*x2
		  //for (int j = 0; j < rasm; j++) //тут я странно сделала много циклов, зочем
		  //	M2[j] = M2[j] * nu0 * sigma;
		  for (int j = 0; j < rasm; j++)
			 b[obj_elems[i].num[j]] += b_loc[j] + M2[j] * nu0 * sigma;
	   }

	   Get_KR1_3D();
	  // z.clear();
	  // L.clear();
	  // D.clear();
	  // U.clear();
	  // z.resize(N);
	  // L.resize(ia[N]);
	  // D.resize(N);
	  // U.resize(ia[N]);
	  // //LOS_LU_3D(t);
   ////	diag_preconditioning(x0_3D[t]);
	  // FactLU(L, U, D);
	  // Direct(L, D, z, b);
	  // Reverse(U, x0_3D[t], z);
	   LOS_LU_3D(t);
	   //diag_preconditioning(x0_3D[1]);

	   out << "t = " << time[t] << endl;
	   //for (int i = 0; i < N; i++) //выводим результат в ребрах
	   //	out << x0_3D[1][i] << endl;
	   //out << endl;

	   f << "t = " << time[t] << endl; //выводим результат в центрах элементов для среза по y

	   for (int i = 0; i < num_elems; i++)
		  if (centers[i].y == centers[center_y].y) //срез по y
		  {
			 otv = Get_solution_3D(t, centers[i].x, centers[i].y, centers[i].z, i, false);
			 otvet.push_back(otv);
			 if (centers[i].x > -30 && centers[i].x < 60 && centers[i].z < 20 && centers[i].z > -60)
				out << centers[i].x << "\t" << centers[i].z << "\t" << otv << endl;
		  }


	   f << endl;

	   //for (int i = 0; i < otvet.size(); i++) //выводим результат в ребрах
	   //   out << centers[i].x << "\t" << centers[i].z << "\t" << otvet[i] << endl;
	   //out << endl;

	   //for (int i = 1; i <= otvet.size() / (num_x - 1); i++) // вывод в файл среза по y
	   //{
		  //for (int j = otvet.size() - (num_x - 1) * i; j < otvet.size() - (num_x - 1) * (i - 1); j++)
			 //f << otvet[j] << "\t";
		  //f << endl;
	   //}
	   //f << endl;

	   //otvet.clear();
	   //finded = false;
	   //out << "t = " << time[t] << endl;
	   //for (int i = 0; i < num_elems && !finded; i++)
	   //{
		  // if (20 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			 //  20 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			 //  -20 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		  // {
			 //  finded = true;
			 //  reciever_elem = i;
		  // }
	   //}

	   //otv = Get_solution_3D(t, 20, 20, -20, reciever_elem, false);
	   //out << "20 20 -20" << "\t" << otv << endl;
	   //finded = false;
	   //for (int i = 0; i < num_elems && !finded; i++)
	   //{
		  // if (0 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			 //  0 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			 //  0 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		  // {
			 //  finded = true;
			 //  reciever_elem = i;
		  // }
	   //}

	   ////f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
	   ////ofstream xout("result_xz.txt");
	   //otv = Get_solution_3D(t, 0, 0, 0, reciever_elem, false);
	   //out << "0 0 0" << "\t" << otv << endl;
	   //finded = false;
	   //for (int i = 0; i < num_elems && !finded; i++)
	   //{
		  // if (50 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			 //  50 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			 //  -50 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		  // {
			 //  finded = true;
			 //  reciever_elem = i;
		  // }
	   //}

	   //f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
	   //ofstream xout("result_xz.txt");
	 /*  otv = Get_solution_3D(t, 50, 50, -50, reciever_elem, false);
	   out << "50 50 -50" << "\t" << otv << endl;*/

	   //ищем значения в приемнике
	    finded = false;
	//   int reciever_elem = 0;

	   for (int i = 0; i < num_elems && !finded; i++)
	   {
		  if (reciever_x >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
			 reciever_y >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
			 reciever_z >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		  {
			 finded = true;
			 reciever_elem = i;
		  }
	   }

	 //  double r = sqrt(reciever_x * reciever_x + reciever_y * reciever_y);
	   //for (int i = 0; i < 4; i++) //прибавляем 2D поле
	   //	x0_3D[1][obj_elems[reciever_elem].num[i]] += Get_solution(1, r, reciever_z, -reciever_y);
	   //for (int i = 4; i < 8; i++)
	   //	x0_3D[1][obj_elems[reciever_elem].num[i]] += Get_solution(1, r, reciever_z, reciever_x);
	   double B2D = Get_solution(t, r, reciever_z, -reciever_y, true);
	   double B3D = Get_solution_3D(t, reciever_x, reciever_y, reciever_z, reciever_elem, true);
	   fout << setprecision(16) << time[t] << "\t" << B2D << "\t" << B2D + B3D << "\t" << B3D << endl;		//ищем B общего поля в приемнике

	   //for (int i = 0; i < 4; i++)
	   //	x0_3D[1][obj_elems[reciever_elem].num[i]] -= Get_solution(1, r, reciever_z, -reciever_y);
	   //for (int i = 4; i < 8; i++)
	   //	x0_3D[1][obj_elems[reciever_elem].num[i]] -= Get_solution(1, r, reciever_z, reciever_x);

	   //fout << setprecision(16) << time[1] << " А отд = " << Get_solution_3D(1, reciever_x, reciever_y, reciever_z, reciever_elem, false) << endl;
	   //fout << setprecision(16) << time[1] << " B2D = " << Get_solution_3D(1, reciever_x, reciever_y, reciever_z, reciever_elem, true) << endl;		//ищем B отд поля в приемнике
	   cout << "2" << endl;
	   //------------------------------------------------------------------------------------------
	}
    //-------------------Инициализация решения на 3-м и посл. слоях-----------------------------

	for (int t = 2; t < time_number; t++)	
	{
		Clear_and_resize(); //очищаем вектора для каждого временного слоя

		//считаем коэффициенты по времени
		t0 = time[t] - time[t - 1]; //t0
		t1 = time[t] - time[t - 2]; //t
		t2 = time[t - 1] - time[t - 2]; //t1
		nu0 = (t1 + t0) / t1 / t0;
		nu1 = t1 / t2 / t0;
		nu2 = t0 / t2 / t1;
		/*nu0 = 1 / t0;
		nu1 = 1 / t0;
		nu2 = 0;*/
		for (int i = 0; i < num_elems; i++)
		{
			Locals_3D(i, t);
			for (int i = 0; i < rasm; i++)
				for (int j = 0; j < rasm; j++)
					A_loc[i][j] = G_loc[i][j] + sigma*M_loc[i][j] * nu0;

			Assemble_Locals_3D(i); //сборка левой части

			//сборка вектора правой части F = b_loc(j)-nu2*M*x2+nu1*M*x1
			for (int j = 0; j < rasm; j++)
				temp[j] = x0_3D[t-2][obj_elems[i].num[j]];
			multiply(M_loc, temp, M2);//M*x2
			for (int j = 0; j < rasm; j++)
				temp[j] = x0_3D[t - 1][obj_elems[i].num[j]];
			multiply(M_loc, temp, M1);//M*x1
			//for (int j = 0; j < rasm; j++) //ненужный цикл
			//{
			//	M1[j] = M1[j] * nu1 * sigma;
			//	M2[j] = M2[j] * nu2 * sigma;
			//}
			for (int j = 0; j < rasm; j++)
				b[obj_elems[i].num[j]] += b_loc[j] - M2[j] * nu2 * sigma + M1[j] * nu1 * sigma;
		}

     	Get_KR1_3D();
	//	z.clear();
	//	L.clear();
	//	D.clear();
	//	U.clear();
	//	z.resize(N);
	//	L.resize(ia[N]);
	//	D.resize(N);
	//	U.resize(ia[N]);
	LOS_LU_3D(t);
	////	diag_preconditioning(x0_3D[t]);
	//	FactLU(L, U, D);
	//	Direct(L, D, z, b);
	//	Reverse(U, x0_3D[t], z);
		//LOS_LU_3D(t);
	//	diag_preconditioning(x0_3D[t]);

		out << "t = " << time[t] << endl;
		//for (int i = 0; i < N; i++) //выводим результат в ребрах
		//	out << x0_3D[t][i] << endl;
		//out << endl;

		//fout << "t = " << time[t] << endl;
		f << "t = " << time[t] << endl; //выводим результат в центрах элементов для среза по y
		//Get_solution_3D(&fout, &f, t);

		for (int i = 0; i < num_elems; i++)
		{
		//	if (centers[i].y == centers[center_y].y) //срез по y
			if (centers[i].y == centers[center_y].y) //срез по y
			{
				otv = Get_solution_3D(t, centers[i].x, centers[i].y, centers[i].z, i, false);
				otvet.push_back(otv);
				if (centers[i].x > -30 && centers[i].x < 60 && centers[i].z < 20 && centers[i].z > -60)
					out << centers[i].x << "\t" << centers[i].z << "\t" << otv << endl;
			}
		}
		//finded = false;
		//out << "t = " << time[t] << endl;
		//for (int i = 0; i < num_elems && !finded; i++)
		//{
		//	if (20 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
		//		20 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
		//		-20 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		//	{
		//		finded = true;
		//		reciever_elem = i;
		//	}
		//}
		//finded = false;
		//otv = Get_solution_3D(t, 20, 20, -20, reciever_elem, false);
		//out << "20 20 -20" << "\t" << otv << endl;

		//for (int i = 0; i < num_elems && !finded; i++)
		//{
		//	if (0 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
		//		0 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
		//		0 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		//	{
		//		finded = true;
		//		reciever_elem = i;
		//	}
		//}
		//finded = false;
		////f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
		////ofstream xout("result_xz.txt");
		//otv = Get_solution_3D(t, 0, 0, 0, reciever_elem, false);
		//out << "0 0 0" << "\t" << otv << endl;

		//for (int i = 0; i < num_elems && !finded; i++)
		//{
		//	if (50 >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
		//		50 >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
		//		-50 >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
		//	{
		//		finded = true;
		//		reciever_elem = i;
		//	}
		//}

		////f << "t = " << time[0] << endl; //выводим результат в центрах элементов для среза по y
		////ofstream xout("result_xz.txt");
		//otv = Get_solution_3D(t, 50, 50, -50, reciever_elem, false);
		//out << "50 50 -50" << "\t" << otv << endl;

		//ищем значения в приемнике
		finded = false;
		int reciever_elem = 0;

		for (int i = 0; i < num_elems && !finded; i++)
		{
			if (reciever_x >= edges[obj_elems[i].num[0]].beg.x && reciever_x <= edges[obj_elems[i].num[0]].end.x &&
				reciever_y >= edges[obj_elems[i].num[4]].beg.y && reciever_y <= edges[obj_elems[i].num[4]].end.y &&
				reciever_z >= edges[obj_elems[i].num[8]].beg.z && reciever_z <= edges[obj_elems[i].num[8]].end.z)
			{
				finded = true;
				reciever_elem = i;
			}
		}
		//for (int i = 0; i < otvet.size(); i++) //выводим результат в ребрах
		//	out << centers[i].x << "\t" << centers[i].z << "\t"<< otvet[i] << endl;
		//out << endl;

		for (int i = 1; i <= otvet.size() / (num_x - 1); i++) // вывод в файл среза по y
		{
			for (int j = otvet.size() - (num_x - 1) * i; j < otvet.size() - (num_x - 1) * (i - 1); j++)
				f << otvet[j] << "\t";
			f << endl;
		}
		f << endl;

		otvet.clear();

		double B2D = Get_solution(t, r, reciever_z, -reciever_y, true);
		double B3D = Get_solution_3D(t, reciever_x, reciever_y, reciever_z, reciever_elem, true);
		fout << setprecision(16) << time[t] << "\t" << B2D << "\t" << B2D + B3D << "\t" << B3D << endl;			//ищем B общего поля в приемнике
		//прибавляем 2D поле
		/*for (int i = 0; i < 4; i++) 
			x0_3D[t][obj_elems[reciever_elem].num[i]] += Get_solution(t, r, reciever_z, -reciever_y);
		for (int i = 4; i < 8; i++)
			x0_3D[t][obj_elems[reciever_elem].num[i]] += Get_solution(t, r, reciever_z, reciever_x);*/

	//	fout << setprecision(16) << time[t] << "\t" << Get_solution_3D(t, reciever_x, reciever_y, reciever_z, reciever_elem, false) << "\t" << Get_solution_3D(t, reciever_x, reciever_y, reciever_z, reciever_elem, true) << endl;

		//for (int i = 0; i < 4; i++) //вычитаем 2D поле
		//	x0_3D[t][obj_elems[reciever_elem].num[i]] -= Get_solution(t, r, reciever_z, -reciever_y);
		//for (int i = 4; i < 8; i++)
		//	x0_3D[t][obj_elems[reciever_elem].num[i]] -= Get_solution(t, r, reciever_z, reciever_x);

		//fout << setprecision(16) << time[t] << " A = " << Get_solution_3D(t, reciever_x, reciever_y, reciever_z, reciever_elem, false) << endl;		//ищем значения поля в приемнике без другого поля
		//fout << setprecision(16) << time[t] << " B2D = " << Get_solution_3D(t, reciever_x, reciever_y, reciever_z, reciever_elem, true) << endl;		//ищем B 2D поля в приемнике
		cout << "3" << endl;
	}
	fout.close();
	out.close();
}

void KURS::Get_KR1_3D() // учет первых краевых
{
		for (int j = 0; j < KR1_3D.size(); j++)
		{
			int node_id = KR1_3D[j];
			di[node_id] = 1E+15;
			b[node_id] = 0;
		/*	for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
				al[k] = 0;
			for (int k = 0; k < ja.size(); k++)
				if (ja[k] == node_id)
					au[k] = 0;*/
			//int node_id = KR1_3D[i][j];
			//di[node_id] = 1;
			//b[node_id] = 0.;
			//for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
			//	al[k] = 0;
			//for (int k = 0; k < ja.size(); k++)
			//	if (ja[k] == node_id)
			//		au[k] = 0;
		}

}
//----------------------------------------------------------------------

//--------------блок функций и процедур для решателя--------------------
void KURS::ATx(vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < N; i++)
	{
		y[i] = di[i] * x[i];
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			int k = ja[j];
			y[i] += au[j] * x[k];
			y[k] += al[j] * x[i];
		}
	}
}

void KURS::Ax(vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < N; i++)
	{
		y[i] = di[i] * x[i];
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			int k = ja[j];
			y[i] += al[j] * x[k];
			y[k] += au[j] * x[i];
		}
	}
}

double KURS::Norm(vector<double>& x)
{
	double norm = 0;
	for (int i = 0; i < N; i++)
		norm += x[i] * x[i];
	return(sqrt(norm));
}

double KURS::mult(const vector<double>& a, const vector<double>& b)
{
	double res = 0;
	for (int i = 0; i < a.size(); i++)
		res += a[i] * b[i];
	return res;
}

void KURS::LOS_LU(int t)
{
	r.resize(N);
	z.resize(N);
	p.resize(N);
	Ar.resize(N);
	y.resize(N);
	L.resize(ia[N]);
	D.resize(N);
	U.resize(ia[N]);
	normB = Norm(b);
	cout.precision(15);
	FactLU(L, U, D);
	double p_p = 0, p_r = 0, r_r = 0, Ar_p = 0;
	double a = 0, B = 0, eps2 = 1e-10;
	Ax(x0[t], y); // y = A * x0
	// y = B - A * x0
	for (int i = 0; i < N; i++)
		y[i] = b[i] - y[i];
	Direct(L, D, r, y); // r0 = L^(-1) * (B - A *x0)
	Reverse(U, z, r); // z0 = U^(-1) * r0
	Ax(z, y); // y = A * z0
	Direct(L, D, p, y); // p0 = L^(-1) * (A *  z0)
	r_r = mult(r, r);
	normR = sqrt(r_r) / normB;
	for (iter = 1; iter < maxIter + 1 && normR >= eps_time; iter++)
	{
		p_p = mult(p, p);
		p_r = mult(p, r);
		a = p_r / p_p;
		// x(k) = x(k-1) + a(k) * z(k-1)
		// r(k) = r(k-1) - a(k) * p(k-1)
		for (int i = 0; i < N; i++)
		{
			x0[t][i] = x0[t][i] + z[i] * a;
			r[i] = r[i] - p[i] * a;
		}
		Reverse(U, y, r); // y = U^(-1) *  r(k)
		Ax(y, Ar); // Ar = A * U^(-1)  *r(k)
		Direct(L, D, Ar, Ar); // Ar = L^(-1) * A *U ^ (-1)* r(k)
		Ar_p = mult(Ar, p); // (Ar, p)
		B = -(Ar_p / p_p);
		// z(k) = U^(-1) * r(k) + B(k) * z(k-1) 
		// p(k) = L^(-1) * A * U^(-1) * r(k) + B(k)* p(k - 1)
		for (int i = 0; i < N; i++)
		{
			z[i] = y[i] + z[i] * B;
			p[i] = Ar[i] + p[i] * B;
		}
		if (r_r - (r_r - a * a * p_p) < eps2)
			r_r = mult(r, r);
		else
			r_r = r_r - a * a * p_p;
		normR = sqrt(r_r) / normB;
		cout << iter << ". " << normR << endl;
	}
}

void KURS::LOS_LU_3D(int t)
{
	double normr0;
	cout << t << endl;
	r.resize(N);
	z.resize(N);
	p.resize(N);
	Ar.resize(N);
	y.resize(N);
	L.resize(ia[N]);
	D.resize(N);
	U.resize(ia[N]);
	normB = Norm(b);
	cout.precision(15);
	FactLU(L, U, D);
	double p_p = 0, p_r = 0, r_r = 0, Ar_p = 0;
	double a = 0, B = 0, eps2 = 1e-10;
	Ax(x0_3D[t], y); // y = A * x0
	// y = B - A * x0
	for (int i = 0; i < N; i++)
		y[i] = b[i] - y[i];
	Direct(L, D, r, y); // r0 = L^(-1) * (B - A *x0)
	Reverse(U, z, r); // z0 = U^(-1) * r0
	Ax(z, y); // y = A * z0
	Direct(L, D, p, y); // p0 = L^(-1) * (A *  z0)
	r_r = mult(r, r);
	normR = sqrt(r_r) / normB;
	normr0 = normR;
	for (iter = 1; iter < maxIter + 1 && normR/normr0 >= eps; iter++)
	{
		p_p = mult(p, p);
		p_r = mult(p, r);
		a = p_r / p_p;
		// x(k) = x(k-1) + a(k) * z(k-1)
		// r(k) = r(k-1) - a(k) * p(k-1)
		for (int i = 0; i < N; i++)
		{
			x0_3D[t][i] = x0_3D[t][i] + z[i] * a;
			r[i] = r[i] - p[i] * a;
		}
		Reverse(U, y, r); // y = U^(-1) *  r(k)
		Ax(y, Ar); // Ar = A * U^(-1)  *r(k)
		Direct(L, D, Ar, Ar); // Ar = L^(-1) * A *U ^ (-1)* r(k)
		Ar_p = mult(Ar, p); // (Ar, p)
		B = -(Ar_p / p_p);
		// z(k) = U^(-1) * r(k) + B(k) * z(k-1) 
		// p(k) = L^(-1) * A * U^(-1) * r(k) + B(k)* p(k - 1)
		for (int i = 0; i < N; i++)
		{
			z[i] = y[i] + z[i] * B;
			p[i] = Ar[i] + p[i] * B;
		}
		if (r_r - (r_r - a * a * p_p) < eps2)
			r_r = mult(r, r);
		else
			r_r = r_r - a * a * p_p;
		normR = sqrt(r_r) / normB;
		cout << iter << ". " << normR/normr0 << endl;
	}
}

void KURS::FactLU(vector<double>& L, vector<double>& U, vector<double>& D)
{
	L = al;
	U = au;
	D = di;
	double l, u, d;
	for (int k = 0; k < N; k++)
	{
		d = 0;
		int i0 = ia[k], i1 = ia[k + 1];
		int i = i0;
		for (; i0 < i1; i0++)
		{
			l = 0;
			u = 0;
			int j0 = i, j1 = i0;
			for (; j0 < j1; j0++)
			{
				int t0 = ia[ja[i0]], t1 = ia[ja[i0] + 1];
				for (; t0 < t1; t0++)
				{
					if (ja[j0] == ja[t0])
					{
						l += L[j0] * U[t0];
						u += L[t0] * U[j0];
					}
				}
			}
			L[i0] -= l;
			U[i0] -= u;
			U[i0] /= D[ja[i0]];
			d += L[i0] * U[i0];
		}
		D[k] -= d;
	}
}

void KURS::Direct(vector<double>& L, vector<double>& D, vector<double>& y, vector<double>& b) // L*y = B
{
	y = b;
	for (int i = 0; i < N; i++)
	{
		double sum = 0;
		int k0 = ia[i], k1 = ia[i + 1];
		int j;
		for (; k0 < k1; k0++)
		{
			j = ja[k0];
			sum += y[j] * L[k0];
		}
		double buf = y[i] - sum;
		y[i] = buf / D[i];
	}
}

void KURS::Reverse(vector<double>& U, vector<double>& x, vector<double>& y) // U*x = y
{
	x = y;
	for (int i = N - 1; i >= 0; i--)
	{
		int k0 = ia[i], k1 = ia[i + 1];
		int j;
		for (; k0 < k1; k0++)
		{
			j = ja[k0];
			x[j] -= x[i] * U[k0];
		}
	}
}

vector<double> Ld, Ll;

void KURS::llt()
{
	int p = 0, m = 0;
	Ll.resize(ia.back());
	Ld.resize(num_edges);
	for (int i = 0; i < num_edges; i++)
	{
		Ld[i] = 0;

		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			Ll[j] = al[j];

			p = ja[j];
			m = ia[i];

			for (int k = ia[p]; k < ia[p + 1]; k++)
			{
				for (int q = m; q < j; q++)
					if (ja[q] == ja[k])
					{
						Ll[j] -= Ll[q] * Ll[k];
						m = q + 1;

						break;
					}
			}

			Ll[j] /= Ld[ja[j]];
			Ld[i] -= Ll[j] * Ll[j];
		}

		Ld[i] = sqrt(di[i] + Ld[i]);
	}
}

void KURS::vec_diff(const vector<double>& x, const vector<double>& y, vector<double>& res)
{
	for (int i = 0; i < x.size(); i++)
		res[i] = x[i] - y[i];
}

void KURS::solve_auxiliary_system(const vector<double>& f, vector<double>& x)
{
	for (int i = 0; i < num_edges; i++)
	{
		double sum = 0.;

		for (int j = ia[i]; j < ia[i + 1]; j++)
			sum += Ll[j] * x[ja[j]];

		x[i] = (f[i] - sum) / Ld[i];
	}

	for (int i = num_edges - 1; i >= 0; i--)
	{
		x[i] /= Ld[i];

		for (int j = ia[i]; j < ia[i + 1]; j++)
			x[ja[j]] -= Ll[j] * x[i];
	}
}

vector <double> s, a;

double KURS::relative_residual(vector<double>& x)
{
	vector<double> rr_vec;
	rr_vec.resize(num_edges);
	Ax(x, rr_vec);
	vec_diff(b, rr_vec, rr_vec);

	return Norm(rr_vec) / Norm(b);
}

void KURS::llt_preconditioning(vector<double>& x0)
{
	cout << "1";
	int total_iter = 1;
	s.resize(num_edges);
	r.resize(num_edges);
	z.resize(num_edges);
	a.resize(num_edges);
	llt();
	if (x0.size() != r.size())
		x0.resize(r.size());
	// r = A * x0
	Ax(x0, r);

	// r = pr - A * x0 = pr - r
	vec_diff(b, r, r);

	// z = M^(-1) * r
	solve_auxiliary_system(r, z);

	double r_norm = Norm(r), rr = r_norm / Norm(b);

	while (rr >= eps && total_iter < maxIter)
	{
		// s = M^(-1) * r
		solve_auxiliary_system(r, s);

		// az = A * z
		Ax(z, a);

		double scal_m_inv_r = mult(s, r), alpha_k = scal_m_inv_r / mult(a, z);

		for (int i = 0; i < num_edges; i++)
		{
			x0[i] += alpha_k * z[i];
			r[i] -= alpha_k * a[i];
		}

		// s = M^(-1) * r
		solve_auxiliary_system(r, s);

		double beta_k = mult(s, r) / scal_m_inv_r;

		for (int i = 0; i < num_edges; i++)
			z[i] = s[i] + beta_k * z[i];

		r_norm = Norm(r);
		rr = r_norm / Norm(b);

		cout << "iter #" << total_iter << " | rr = " << rr << endl;
		total_iter++;
	}
	cout << "RELATIVE RESIDUAL = " << relative_residual(x0) << endl;
}

void KURS::diag_preconditioning(vector<double>& x0)
{
	s.resize(num_edges);
	r.resize(num_edges);
	z.resize(num_edges);
	az.resize(num_edges);
	// r = A * x0
	Ax(x0, r);

	// r = pr - A * x0 = pr - r
	vec_diff(b, r, r);

	// z = M^(-1) * r
	for (int i = 0; i < num_edges; i++)
		z[i] = r[i] / di[i];
	int total_iter = 1;
	double r_norm = Norm(r), rr = r_norm / Norm(b);

	while (rr >= eps && total_iter < maxIter)
	{
		// az = A * z
		Ax(z, az);

		double scal_m_inv_r = mult(z, r), alpha_k = scal_m_inv_r / mult(az, z);

		for (int i = 0; i < num_edges; i++)
		{
			x0[i] += alpha_k * z[i];
			r[i] -= alpha_k * az[i];
		}

		// s = M^(-1) * r
		for (int i = 0; i < num_edges; i++)
			s[i] = r[i] / di[i];

		double beta_k = mult(s, r) / scal_m_inv_r;

		for (int i = 0; i < num_edges; i++)
			z[i] = s[i] + beta_k * z[i];

		r_norm = Norm(r);
		rr = r_norm / Norm(b);

		cout << "iter #" << total_iter << " | rr = " << rr << endl;

		total_iter++;
	}

	cout << "RELATIVE RESIDUAL = " << relative_residual(x0) << endl;
}


//----------------------------------------------------------------------

int main()
{
	KURS A;
	A.Input();

	A.rasm = 4; //для двумерной
	A.Generate_Portrait();
	A.Find_field();

	A.rasm = 12; //для трехмерной
	A.Generate_Portrait_3D();
	A.Find_in_object();
	return 0;
}
