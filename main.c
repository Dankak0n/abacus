#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <malloc.h>

#define SIZE 400
#define EPS 0.00001
#define INTERNAL_FUNCTIONS_CNT 18

char file_name[] = "input.txt";
enum ftype_t{LEFTRIGHT, RIGHT, LEFT};

typedef struct{
    int l, r;
} pair;

typedef struct{
    float real, imag;
} number_t;

typedef struct{
    char name[SIZE];
    int priority;
    enum ftype_t type;
    int argc;
    number_t (*function_ptr)(int, number_t *);
} function_opa;

function_opa internal_functinos[INTERNAL_FUNCTIONS_CNT];
int known_values_cnt;

int get_priority(char * str) {
    for (int i = 0; i < INTERNAL_FUNCTIONS_CNT; i++) {
        if (!strncmp(str, internal_functinos[i].name, strlen(internal_functinos[i].name))) {
            return internal_functinos[i].priority;
        }
    }
    return -1000;
}

enum ftype_t get_type(char * str) {
    for (int i = 0; i < INTERNAL_FUNCTIONS_CNT; i++) {
        if (!strncmp(str, internal_functinos[i].name, strlen(internal_functinos[i].name))) {
            return internal_functinos[i].type;
        }
    }
}

int get_argc(char * str) {
    for (int i = 0; i < INTERNAL_FUNCTIONS_CNT; i++) {
        if (!strncmp(str, internal_functinos[i].name, strlen(internal_functinos[i].name))) {
            return internal_functinos[i].argc;
        }
    }
}

number_t get_value(char * str, number_t * argv) {
    for (int i = 0; i < INTERNAL_FUNCTIONS_CNT; i++) {
        if (!strncmp(str, internal_functinos[i].name, strlen(internal_functinos[i].name))) {
            return internal_functinos[i].function_ptr(get_argc(str), argv);
        }
    }
}

typedef struct{
    char name[SIZE];
    number_t value;
} defined_value_t;

defined_value_t known_values[SIZE];

number_t ADD(int argc, number_t * argv);

number_t SUB(int argc, number_t * argv);

number_t MUL(int argc, number_t * argv);

number_t DIV(int argc, number_t * argv);

number_t SIN(int argc, number_t * argv);

number_t COS(int argc, number_t * argv);

number_t TAN(int argc, number_t * argv);

number_t LOG(int argc, number_t * argv);

number_t LGE(int argc, number_t * argv);

number_t SQR(int argc, number_t * argv);

number_t POW(int argc, number_t * argv);

number_t ABS(int argc, number_t * argv);

number_t EXP(int argc, number_t * argv);

number_t REL(int argc, number_t * argv);

number_t IMG(int argc, number_t * argv);

number_t MAG(int argc, number_t * argv);

number_t PHS(int argc, number_t * argv);

number_t IMJ(int argc, number_t * argv);

int is_reserved(char symb);

int is_to_number(char symb);

void push(pair ** arr, int * size, pair x);

int parse(char (* str)[SIZE], pair ** arr);

void init();

number_t get_simplest(char * str, int size);

void get_argv(number_t * argv, char * str, pair * arr, int l, int r, int pos_max, int argc, enum ftype_t type);

number_t go_(char * str, pair * arr, int l, int r);

int main() {
    init();
#ifdef LOCAL
    freopen(file_name, "r", stdin);
#endif
    char expression_str[SIZE][SIZE];
    pair * expression_arr[SIZE];
    int size[SIZE], id = 0;
    while (gets(expression_str[id])) {
        size[id] = parse(&expression_str[id], &expression_arr[id]);
        id++;
    }
    while (--id) {
        strncpy(known_values[known_values_cnt].name, expression_str[id], expression_arr[id][0].r);
        known_values[known_values_cnt].value = go_(expression_str[id], expression_arr[id], 2, size[id]);
        known_values_cnt++;
    }
    number_t ans = go_(expression_str[0], expression_arr[0], 0, size[0]);
    printf("%f + %fj", ans.real, ans.imag);
}

void init() {
    internal_functinos[0]  = (function_opa){"+",     200, LEFTRIGHT, 2, ADD};
    internal_functinos[1]  = (function_opa){"-",     200, LEFTRIGHT, 2, SUB};
    internal_functinos[2]  = (function_opa){"*",     100, LEFTRIGHT, 2, MUL};
    internal_functinos[3]  = (function_opa){"/",     100, LEFTRIGHT, 2, DIV};
    internal_functinos[4]  = (function_opa){"sin",   0,   RIGHT,     1, SIN};
    internal_functinos[5]  = (function_opa){"cos",   0,   RIGHT,     1, COS};
    internal_functinos[6]  = (function_opa){"tg",    0,   RIGHT,     1, TAN};
    internal_functinos[6]  = (function_opa){"log",   0,   RIGHT,     2, LOG};
    internal_functinos[7]  = (function_opa){"ln",    0,   RIGHT,     1, LGE};
    internal_functinos[8]  = (function_opa){"sqrt",  0,   RIGHT,     1, SQR};
    internal_functinos[9]  = (function_opa){"pow",   0,   RIGHT,     2, POW};
    internal_functinos[10] = (function_opa){"abs",   0,   RIGHT,     1, ABS};
    internal_functinos[11] = (function_opa){"exp",   0,   RIGHT,     1, EXP};
    internal_functinos[12] = (function_opa){"real",  0,   RIGHT,     1, REL};
    internal_functinos[13] = (function_opa){"imag",  0,   RIGHT,     1, IMG};
    internal_functinos[14] = (function_opa){"mag",   0,   RIGHT,     1, MAG};
    internal_functinos[15] = (function_opa){"phase", 0,   RIGHT,     1, PHS};
    internal_functinos[16] = (function_opa){"^",     50,  LEFTRIGHT, 2, POW};
    internal_functinos[17] = (function_opa){"j",     0,   LEFT,      1, IMJ};

    known_values[0]  = (defined_value_t){"PI", (number_t){M_PI, 0}};
    known_values[1]  = (defined_value_t){"e",  (number_t){M_E, 0}};
    known_values_cnt = 2;
}

int parse(char (* str)[SIZE], pair ** arr) {
    char kukozhed[SIZE];
    pair * ans = NULL;
    int size = 0, ans_size = 0;
    for (int i = 0; (*str)[i] != '\0'; i++) {
        if ((*str)[i] != ' ' && (*str)[i] != '\n') {
            kukozhed[size] = (*str)[i];
            size++;
        }
    }
    kukozhed[size] = '\0';
    strcpy(*str, kukozhed);
    for (int i = 0; i < size; i++) {
        if (is_reserved(kukozhed[i])) {
            push(&ans, &ans_size, (pair){i, i + 1});
        } else if (is_to_number(kukozhed[i])) {
            int l = i;
            i++;
            while (i < size && is_to_number(kukozhed[i])) i++;
            push(&ans, &ans_size, (pair){l, i});
            i--;
        } else {
            int l = i;
            i++;
            while (i < size && (!is_reserved(kukozhed[i]) || kukozhed[i] == 'j')) i++;
            push(&ans, &ans_size, (pair){l, i});
            i--;
        }
    }
    *arr = ans;
    return ans_size;
}

number_t get_simplest(char * str, int size) {
    for (int i = 0; i < known_values_cnt; i++) {
        if (!strncmp(str, known_values[i].name, strlen(known_values[i].name))) {
            return known_values[i].value;
        }
    }
    number_t ans = (number_t){0, 0};
    int was_dot = 0;
    for (int i = 0; i < size; i++) {
        if (str[i] == '.') {
            was_dot = 1;
            continue;
        }
        if (!was_dot) {
            ans.real *= 10.0;
            ans.real += str[i] - '0';
        }
        if (was_dot) {
            ans.real += (float)(str[i] - '0') * pow(10.0, (float)-was_dot);
            was_dot++;
        }
    }
    return ans;
}

void get_argv(number_t * argv, char * str, pair * arr, int l, int r, int pos_max, int argc, enum ftype_t type) {
    if (type == LEFTRIGHT) {
        argv[0] = go_(str, arr, l, pos_max);
        argv[1] = go_(str, arr, pos_max + 1, r);
    }
    if (type == RIGHT) {
        l = pos_max + 2;
        int nesting = 0, left = l, size = 0;
        for (int i = l; -1 <= nesting; i++) {
            if (str[arr[i].l] == '(' || str[arr[i].l] == ')') {
                if (str[arr[i].l] == '(')
                    nesting++;
                else
                    nesting--;
            }
            if (nesting < 1) {
                if (str[arr[i].l] == ',' || (str[arr[i].l] == ')' && nesting == -1)) {
                    argv[size] = go_(str, arr, left, i);
                    left = i + 1;
                    size++;
                }
            }
            if (nesting == -1) break;
        }
    }
    if (type == LEFT) {
        if (r - l == 1) {
            argv[0] = (number_t){1.0, 0};
        } else {
            argv[0] = (number_t){go_(str, arr, l, r - 1).real, 0};
        }
    }
}

number_t go_(char * str, pair * arr, int l, int r) {
    if (r - l < 1) {
        return (number_t){0, 0};
    }
    if (r - l == 1 && str[arr[l].l] != 'j') {
        return get_simplest(str + arr[l].l, arr[l].r - arr[l].l);
    }
    int nesting = 0, max_priority = -100, pos_max = -1, cnt_zero = 0;
    for (int i = l; i < r; i++) {
        if (str[arr[i].l] == '(' || str[arr[i].l] == ')') {
            if (str[arr[i].l] == '(')
                nesting++;
            else
                nesting--;
            if (nesting == 0) cnt_zero++;
            continue;
        }
        if (nesting == 0) {
            if (get_priority(str + arr[i].l) >= max_priority) {
                max_priority = get_priority(str + arr[i].l);
                pos_max = i;
            }
        }
    }
    if (cnt_zero == 1 && str[arr[l].l] == '(' && str[arr[r - 1].l] == ')') {
        return go_(str, arr, l + 1, r - 1);
    }
    char * last_func = str + arr[pos_max].l;
    enum ftype_t type = get_type(last_func);
    int argc = get_argc(last_func);
    number_t * argv = (number_t *)malloc(sizeof(number_t) * argc);
    get_argv(argv, str, arr, l, r, pos_max, argc, type);
    return get_value(last_func, argv);
}

number_t ADD(int argc, number_t * argv) {
    return (number_t){
        argv[0].real + argv[1].real,
        argv[0].imag + argv[1].imag
    };
}

number_t SUB(int argc, number_t * argv) {
    return (number_t){
        argv[0].real - argv[1].real,
        argv[0].imag - argv[1].imag
    };
}

number_t MUL(int argc, number_t * argv) {
    return (number_t){
        argv[0].real * argv[1].real - argv[0].imag * argv[1].imag,
        argv[0].imag * argv[1].real + argv[0].real * argv[1].imag
    };
}

number_t DIV(int argc, number_t * argv) {
    return (number_t){
        (argv[0].real * argv[1].real + argv[0].imag * argv[1].imag) /
        (pow(argv[1].real, 2) + pow(argv[1].imag, 2)),
        (argv[0].imag * argv[1].real - argv[0].real * argv[1].imag) /
        (pow(argv[1].real, 2) + pow(argv[1].imag, 2)),
    };
}

number_t SIN(int argc, number_t * argv) {
    return (number_t){sin(argv[0].real), 0};
}

number_t COS(int argc, number_t * argv) {
    return (number_t){cos(argv[0].real), 0};
}

number_t TAN(int argc, number_t * argv) {
    return (number_t){tan(argv[0].real), 0};
}

number_t LOG(int argc, number_t * argv) {
    return (number_t){log(argv[0].real) / log(argv[1].real), 0};
}

number_t LGE(int argc, number_t * argv) {
    return (number_t){log(argv[0].real), 0};
}

number_t SQR(int argc, number_t * argv) {
    if (fabs(argv[0].real) < EPS && fabs(argv[0].imag) < EPS) {
        return (number_t){0, 0};
    }
    number_t phase, mod;
    if (argv[0].real < -EPS && fabs(argv[0].imag) < EPS) {
        phase = (number_t){M_PI, 0};
    } else {
        phase = (number_t){2 * atan(argv[0].imag / (sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)) + argv[0].real)), 0};
    }
    mod = (number_t){sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)), 0};
    return (number_t){pow(mod.real, 0.5) * cos(0.5 * phase.real), pow(mod.real, 0.5) * sin(0.5 * phase.real)};
}

number_t POW(int argc, number_t * argv) {
    if (fabs(argv[0].real) > EPS) {
        number_t phase, mod;
        if (argv[0].real < -EPS && fabs(argv[0].imag) < EPS) {
            phase = (number_t){M_PI, 0};
        } else {
            phase = (number_t){2 * atan(argv[0].imag / (sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)) + argv[0].real)), 0};
        }
        mod = (number_t){sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)), 0};
        return (number_t){pow(mod.real, argv[1].real) * cos(argv[1].real * phase.real), pow(mod.real, argv[1].real) * sin(argv[1].real * phase.real)};
    } else {
        return (number_t){pow(argv[0].real, argv[1].real) * cos(log(argv[0].real) * argv[1].imag), pow(argv[0].real, argv[1].real) * sin(log(argv[0].real) * argv[1].imag)};
    }
}

number_t ABS(int argc, number_t * argv) {
    return (number_t){sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)), 0};
}

number_t EXP(int argc, number_t * argv) {
    return (number_t){exp(argv[0].real) * cos(argv[0].imag), exp(argv[0].real) * sin(argv[0].imag)};;
}

number_t REL(int argc, number_t * argv) {
    return (number_t){argv[0].real, 0};
}

number_t IMG(int argc, number_t * argv) {
    return (number_t){argv[0].imag, 0};
}

number_t MAG(int argc, number_t * argv) {
    return (number_t){sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)), 0};
}

number_t PHS(int argc, number_t * argv) {
    if (argv[0].real < -EPS && fabs(argv[0].imag) < EPS) {
        return (number_t){M_PI, 0};
    } else {
        return (number_t){2 * atan(argv[0].imag / (sqrt(pow(argv[0].real, 2) + pow(argv[0].imag, 2)) + argv[0].real)), 0};
    }
}

number_t IMJ(int argc, number_t * argv) {
    return (number_t){0, argv[0].real};
}

char * reserved = "()+-*/^j,=";

int is_reserved(char symb) {
    for (int i = 0; reserved[i] != '\0'; i++) {
        if (reserved[i] == symb)
            return 1;
    }
    return 0;
}

int is_to_number(char symb) {
    if ('0' <= symb && symb <= '9') return 1;
    if (symb == '.') return 1;
    return 0;
}

void push(pair ** arr, int * size, pair x) {
    (*size)++;
    pair * tmp = (pair *)malloc((*size) * sizeof(pair));
    for (int i = 0; i < (*size); i++) {
        tmp[i] = (i < (*size) - 1 ? (*arr)[i] : x);
    }
    *arr = tmp;
}

