#include <bits/stdc++.h>
using namespace std;

#define PI acos(-1.0)
#define delta 0.98
#define ITERS 10
#define T_init 10
#define T_min 1
#define INF -1e9 

double k1[ITERS], k2[ITERS], k3[ITERS];
int F_old[ITERS];

const int N = 1000;
const int max_guesses = 1e5;
int guess_num = 0;

unordered_map<string, int> LDS;
unordered_map<string, int> L;
unordered_map<string, int> D;
unordered_map<string, int> S;
unordered_map<int, int> cntLDS;
unordered_map<int, int> cntL;
unordered_map<int, int> cntD;
unordered_map<int, int> cntS;
unordered_set<string> Test;
unordered_set<string> Train;
unordered_set<string> cracked_all;

unordered_map<string, int> mapCracked;
int produced;
int cracked;
ofstream ofs_XY;

struct pNode
{
    string s;
    double p;
    bool operator< (const pNode &t) const
    {
        return p - t.p > 1e-10;
    }
};

struct Node
{
    string s;
    string lds;
    int pivot;
    double p;
    vector<int> idx;
};
queue<Node*> guessQue;      // 预处理得到的口令集

struct cmp
{
    bool operator () (Node* &a, Node* &b) const
    {
        return a->p - b->p < 1e-10;
    }
};

// 存储每个字符串的概率
vector<pNode> vecL[N];
vector<pNode> vecD[N];
vector<pNode> vecS[N];

static int getNum(string &s, int i, int len)
{
    int num = 0;
    while (i < len && isdigit(s[i]))
    {
        num = num*10 + s[i] - '0';
        ++i;
    }
    return num;
}

static Node* getFirstNode(string lds, double c1, double c2, double c3)
{
    double p = LDS[lds] * 1.0 / cntLDS[lds.size()];
    int len = lds.size(), num;
    vector<int> idx;
    string s;

    for (int i = 0; i < len; i++)
    {
        if (isdigit(lds[i])) continue;
        num = getNum(lds, i+1, len);
        if (lds[i] == 'L')
        {
            idx.push_back(1);
            s += vecL[num][0].s;
            p *= vecL[num][0].p * c1;
        }
        else if (lds[i] == 'D')
        {
            idx.push_back(1);
            s += vecD[num][0].s;
            p *= vecD[num][0].p * c2;
        }
        else
        {
            idx.push_back(1);
            s += vecS[num][0].s;
            p *= vecS[num][0].p * c3;
        }
    }

    Node* tmp = new Node;
    tmp->idx = idx;
    tmp->p = p;
    tmp->s = s;
    tmp->lds = lds;
    tmp->pivot = 0;
    return tmp;
}

static void getNextNode(priority_queue<Node*, vector<Node*>, cmp> &tmpQue, Node* f, double c1, double c2, double c3)
{
    int k = 0, kk, idx, num;
    for (int i = 0; i < f->pivot; ++i)
    {
        while (isdigit((f->lds)[++k]));
    }
    k--;
    
    string s;
    char c;
    for (int i = f->pivot; i < (int)(f->idx).size(); ++i)
    {
        if(guessQue.size() + tmpQue.size() >= max_guesses) return ;

        idx = (f->idx)[i];
        while (isdigit((f->lds)[++k]));
        num = getNum(f->lds, k+1, (f->lds).size());
        c = (f->lds)[k];
        if (c == 'L' && (int)vecL[num].size() == idx) continue;
        else if (c == 'D' && (int)vecD[num].size() == idx) continue;
        else if (c == 'S' && (int)vecS[num].size() == idx) continue;

        Node *tmp = new Node;
        tmp->lds = f->lds;
        tmp->p = f->p;
        if (c == 'L')
        {
            tmp->p /= vecL[num][idx-1].p;
            tmp->p *= vecL[num][idx].p * c1;
        }
        else if (c == 'D')
        {
            tmp->p /= vecD[num][idx-1].p;
            tmp->p *= vecD[num][idx].p * c2;
        }
        else
        {
            tmp->p /= vecS[num][idx-1].p;
            tmp->p *= vecS[num][idx].p * c3;
        }
        tmp->idx = f->idx;
        tmp->idx[i]++;

        s.clear();
        kk = -1;
        for (int j = 0; j < (int)(f->idx).size(); j++)
        {
            while (isdigit((f->lds)[++kk]));
            num = getNum(f->lds, kk+1, (f->lds).size());
            if ((tmp->lds)[kk] == 'L') s += vecL[num][(tmp->idx)[j]-1].s;
            else if ((tmp->lds)[kk] == 'D') s += vecD[num][(tmp->idx)[j]-1].s;
            else s += vecS[num][(tmp->idx)[j]-1].s;
        }
        tmp->s = s;
        tmp->pivot = i;
        tmpQue.push(tmp);
    }
}

// 读取测试集
static void getTestData(int option)
{
    ifstream ifs;
    if(option == 1) ifs.open("yahoo_test1.txt", ios::in);
    else if(option == 2) ifs.open("yahoo_test2.txt", ios::in);
    if (!ifs.is_open()) cout << "yahoo_test.txt open fail!" << endl;

    // 读取所有的测试数据
    string buf;
    Test.clear();
    while (getline(ifs, buf)) Test.insert(buf);
    ifs.close();
}

// 测试
static int test(double c1, double c2, double c3, int option)
{
    priority_queue<Node*, vector<Node*>, cmp> testQue;
    queue<Node*> tmpQue;
    produced = 0;
    cracked = 0;

    while(!guessQue.empty())
    {
        Node *node = new Node;
        Node *tmp = guessQue.front();
        node->idx = tmp->idx;
        node->lds = tmp->lds;
        node->p = tmp->p;
        node->pivot = tmp->pivot;
        node->s = tmp->s;
        guessQue.pop();

        int len = (node->lds).length();
        for(int i = 0; i < len; ++i)
        {
            if(isdigit((node->lds)[i])) continue;

            if((node->lds)[i] == 'L') node->p *= c1;
            else if((node->lds)[i] == 'D') node->p *= c2;
            else if((node->lds)[i] == 'S') node->p *= c3;
        }

        testQue.push(node);
        tmpQue.push(tmp);
    }
    while(!tmpQue.empty())
    {
        Node *tmp = tmpQue.front();
        tmpQue.pop();
        guessQue.push(tmp);
    }

    mapCracked.clear();
    while (!testQue.empty()) 
    {
        Node *f = testQue.top();
        testQue.pop();
        produced++;

        if (Test.count(f->s) && !mapCracked[f->s])
        {
            cracked++;
            mapCracked[f->s] = 1;
        }
        delete f;       // 删除节点，释放空间

        if(produced >= guess_num)
        {
            while(!testQue.empty())
            {
                f = testQue.top();
                testQue.pop();
                delete f; 
            }
            break;
        }
    }
    
    if(option == 1) ofs_XY.open("orignal.txt", ios::app | ios :: out);
    else if(option == 2) ofs_XY.open("optimal.txt", ios::app | ios :: out);
    if(option) 
    {
        ofs_XY << produced << "," <<  cracked << endl;
        ofs_XY.close();
    }
    return cracked;
}

// 预先生成大量猜测集
static void genGuessSet()
{
    priority_queue<Node*, vector<Node*>, cmp> tmpQue;
    double c1 = 1.0, c2 = 1.0, c3 = 1.0;

    for (auto it : Train)
    {
        if (LDS.find(it) == LDS.end() || tmpQue.size() >= max_guesses) continue;
        Node *tmp = getFirstNode(it, c1, c2, c3);
        tmpQue.push(tmp);
    }

    while(!tmpQue.empty())
    {
        Node *f = tmpQue.top();
        tmpQue.pop();
        guessQue.push(f);

        if(guessQue.size() + tmpQue.size() < max_guesses) getNextNode(tmpQue, f, c1, c2, c3);
    }
    printf("%d\n", guessQue.size());
}

// 初始化每个概率
static void initVec()
{
    ifstream ifs;
    pNode tmp;
    string buf;
    double p;
    int len;
    int num;

    ifs.open("XZB/yahoo_L_INFLUENCE.txt", ios::in);
    if (!ifs.is_open()) cout << "XZB_yahoo_train_output_L.txt open fail!" << endl;
    while (ifs >> buf >> p)
    {
        len = buf.size();
        tmp = {buf, p};
        vecL[len].push_back(tmp);
    }
    ifs.close();

    ifs.open("XZB/yahoo_D_INFLUENCE.txt", ios::in);
    if (!ifs.is_open()) cout << "XZB_yahoo_train_output_D.txt open fail!" << endl;
    while (ifs >> buf >> p) 
    {
        len = buf.size();
        tmp = {buf, p};
        vecD[len].push_back(tmp);
    }
    ifs.close();

    ifs.open("XZB/yahoo_S_INFLUENCE.txt", ios::in);
    if (!ifs.is_open()) cout << "XZB_yahoo_train_output_S.txt open fail!" << endl;
    while (ifs >> buf >> p)
    {
        len = buf.size();
        tmp = {buf, p};
        vecS[len].push_back(tmp);
    }
    ifs.close();

    for (int i = 0; i < N; ++i)
    {
        if (vecL[i].empty()) continue;
        sort(vecL[i].begin(), vecL[i].end());
    }
    for (int i = 0; i < N; ++i)
    {
        if (vecD[i].empty()) continue;
        sort(vecD[i].begin(), vecD[i].end());
    }
    for (int i = 0; i < N; ++i)
    {
        if (vecS[i].empty()) continue;
        sort(vecS[i].begin(), vecS[i].end());
    }
}

// 获得训练集
static void train()
{
    ifstream ifs;
    ifs.open("LDS_yahoo_train.txt", ios::in);
    if (!ifs.is_open()) cout << "LDS_yahoo_train.txt open fail!" << endl;

    string buf;
    int num;
    while (ifs >> buf >> num)
    {
        LDS[buf] = num;
        cntLDS[(int)buf.size()] += num;
        Train.insert(buf);  
    }

    ifs.close();
    initVec();
    cout << "training has done!" << endl;
}

static void init()
{
    srand((unsigned)time(NULL));
    k1[0] = 1.0, k2[1] = 1.0, k3[0] = 1.0;
	for(int i = 1; i < ITERS; ++i)
		k1[i] = (double(rand() % 1000001)) / 1000000.0 * 3.0;
	for(int i = 1; i < ITERS; ++i)
		k2[i] = (double(rand() % 1000001)) / 1000000.0 * 3.0;
    for(int i = 1; i < ITERS; ++i)
		k3[i] = (double(rand() % 1000001)) / 1000000.0 * 3.0;
}

static void SA()
{
    ofstream ofs_k;
    srand((unsigned)time(NULL));
	double T = T_init;                          // 温度初始值 
	int F_max = test(1.0, 1.0, 1.0, 0);         // 结果初始为最小 
	double c1 = 1.0, c2 = 1.0, c3 = 1.0;
    printf("original_PCFG: %d\n", F_max);

	while(T > T_min) 
    {
		for(int i = 0; i < ITERS; ++i) 
        { 
			double new_k1 = k1[i] + (double(rand() % 1001)) / 100000.0 * T * (rand() % 2 ? -1.0 : 1.0);
			double new_k2 = k2[i] + (double(rand() % 1001)) / 100000.0 * T * (rand() % 2 ? -1.0 : 1.0);
            double new_k3 = k3[i] + (double(rand() % 1001)) / 100000.0 * T * (rand() % 2 ? -1.0 : 1.0);

			if(new_k1 < 0.001 || new_k1 > 3.0) continue; 
			if(new_k2 < 0.001 || new_k2 > 3.0) continue;
			if(new_k3 < 0.001 || new_k3 > 3.0) continue;

			int F_new = test(new_k1, new_k2, new_k3, 0);
			if(F_new > F_old[i]) 
            {  
                F_old[i] = F_new;
                k1[i] = new_k1;
				k2[i] = new_k2;
                k3[i] = new_k3;

                if(F_new > F_max)
                {
                    F_max = F_new;
                    c1 = new_k1;
                    c2 = new_k2;
                    c3 = new_k3;
                }
			}
            else
            { 
				double P = exp((double)(F_new - F_old[i]) / (T * 10.0));  // 最大值/最小值问题，计算该概率 P 公式不同 
				if(P > (double(rand() % 101)) / 100.0) 
                {
					F_old[i] = F_new;
                    k1[i] = new_k1;
                    k2[i] = new_k2;
                    k3[i] = new_k3;
				}
			}
		}
		T = T * delta; 
	}

    ofs_k.open("constant.txt", ios::out | ios::app);
    ofs_k << c1 << " " << c2 << " " << c3 << endl;
    ofs_k.close();
    printf("optimal_PCFG: %d\n", F_max);
}

static void result_test()
{
    ifstream ifs;
    int index = 0;
    guess_num = 0;
    double c1[N], c2[N], c3[N];

    getTestData(2);
    ifs.open("constant.txt", ios::in);
    while(ifs >> c1[++index] >> c2[index] >> c3[index]);
    ifs.close();

    printf("index : %d\n", index);
    for(int i = 1; i < index; ++i)
    {
        guess_num += 100;
        test(1.0, 1.0, 1.0, 1);
        test(c1[i], c2[i], c3[i], 2);
    }
}

int main()
{
    srand((unsigned)time(NULL));
    ios::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);

    train();
    getTestData(1);
    genGuessSet();

    while(guess_num < 2000) 
    {
        guess_num += 100;
        init();
        SA();
    }
    
    result_test();
    return 0;
}
