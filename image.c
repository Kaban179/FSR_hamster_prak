#include "lodepng.c"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define eps 4
#define pi 3.14159265




char* load_png_file(const char *filename, int *width, int *height) {
    unsigned char *image = NULL;
    int error = lodepng_decode32_file(&image, width, height, filename);
    if (error) {
        printf("error %u: %s\n", error, lodepng_error_text(error));
        return NULL;
    }

    return (image);
}


typedef struct Node {
    int size;
    unsigned char R, G, B, alpha;
    struct Node *up, *down, *left, *right;
    struct Node *par;

}Node;


unsigned char max(unsigned char p1, unsigned char p2){
    if (p1 > p2){
        return p1;
    }
    return p2;
}


unsigned char min(unsigned char p1, unsigned char p2){
    if (p1 < p2){
        return p1;
    }
    return p2;
}


void make_picture_grey(unsigned char* image, int w, int h){
    unsigned char R, G, B, alpha;
    for (int i = 0; i < h * w; i++){
        R = image[i * 4];
        G = image[i * 4 + 1];
        B = image[i * 4 + 2];
        alpha = image[i * 4 + 3];
        unsigned char grey = (R + G + B) / 3;
        image[4 * i] = grey;
        image[4 * i + 1] = grey;
        image[4 * i + 2] = grey;
        image[4 * i + 3] = 255;
    }
    return;
}


void color_pixel(unsigned char* image, int w, int h, int i, unsigned char R, unsigned char G, unsigned char B){
    image[4 * i] = R;
    image[4 * i + 1] = G;
    image[4 * i + 2] = B;
    image[4 * i + 3] = 255;
    return;
}


void make_Gauss_matrix(double matrix[5][5], double Sigma){

    int i, j;
    for (i = -2; i < 3; i++){
        for (j = -2; j < 3; j++){
            matrix[i + 2][j + 2] = (1 / (2 * pi * Sigma * Sigma)) * (exp(-(i * i + j * j) / (2 * Sigma * Sigma)));
        }
    }
    return;
}


void Gauss_filter(unsigned char* picture, int w, int h, double Sigma){
    int i, j, dx, dy;
    unsigned char* new_image = malloc(w * h * 4 * sizeof(unsigned char));
    double Gauss_matrix [5][5];
    double R = 0.0;
    make_Gauss_matrix(Gauss_matrix, Sigma);
    printf("Gauss matrix with sigma = %lf\n", Sigma);
    for (i = -2; i < 3; i++){
        for (j = -2; j < 3; j++){
            printf("%lf ", Gauss_matrix[i + 2][j + 2]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 2; i < w - 2; i++){
        for (j = 2; j < h - 2; j++){
            R = 0.0;
            for (dx = -2; dx < 3; dx++){
                for (dy = -2; dy < 3; dy++){
                    R += (Gauss_matrix[dx + 2][dy + 2]) * (picture[((j+dy) * w + (i+dx)) * 4]);
                }
            }
            new_image[(j * w + i) * 4] = (int)round(R);
            new_image[(j * w + i) * 4 + 1] =  (int)round(R);
            new_image[(j * w + i) * 4 + 2] =  (int)round(R);
            new_image[(j * w + i) * 4 + 3] = 255;
        }
    }
    char* filename = "image_after_Gauss.png";
    for (i = 0; i < w * h * 4; i++){
        picture[i] = new_image[i];
    }
    lodepng_encode32_file(filename, new_image, w, h);
    free(new_image);
}


void applySobelFilter(unsigned char *image, int w, int h){
    int gx[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
    int gy[3][3] = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
    unsigned char *result = malloc(w * h * 4 * sizeof(unsigned char));
    unsigned char *result2 = malloc(w * h * 4 * sizeof(unsigned char));
    for (int y = 1; y < h - 1; y++) {
        for (int x = 1; x < w - 1; x++) {
            int sumX = 0, sumY = 0;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int index = ((y+dy) * w + (x+dx)) * 4;
                    int gray = (image[index] + image[index+1] + image [index+2]) / 3;
                    sumX += gx[dy+1][dx+1] * gray;
                    sumY += gy[dy+1][dx+1] * gray;
                }
            }
            int magnitude = sqrt(sumX * sumX + sumY * sumY);
            if (magnitude > 255) magnitude = 255;
            if (magnitude < 0) magnitude = 0;
            int resultIndex = (y * w + x) * 4;
            result[resultIndex] = (unsigned char) magnitude;
            result[resultIndex+1] = (unsigned char) magnitude;
            result[resultIndex+2] = (unsigned char) magnitude;
            result[resultIndex+3] = image[resultIndex+3];
        }
    }
    for (int i = 0; i < w * h * 4; i++) {
        result2[i] = result[i];
    }
    char* filename = "image_after_Sobel.png";
    for (int i = 0; i < w * h * 4; i++) {
        result[i] = result2[i] * 0.8 + image[i] * 0.7;
        if (result2[i] * 0.8 + image[i] * 0.7 > 255){
            result[i] = 255;
        }
    }
    for (int i = 0; i < w * h * 4; i++) {
        if (result2[i] * 0.8 + image[i] * 0.7 > 255){
            image[i] = 255;
        }
        else{
            image[i] = result2[i] * 0.8 + image[i] * 0.7;
        }
    }
    lodepng_encode32_file(filename, result, w, h);
    free(result);
    free(result2);
}



Node* find_set(Node* el) {
    while (el -> par != el){
        el = el -> par;
    }
    return el -> par;
}

void union_set(Node* node1, Node* node2) {
    if (node1 && node2){
        Node* par1 = find_set(node1);
        Node* par2 = find_set(node2);
        double color_difference = fabs(node1 -> R - node2 -> R);
        if (node1 -> R < 10) return;
        if (par1 != par2){
            if (color_difference < eps){
                if (par1 -> size > par2 -> size){
                    par2 -> par = par1;
                }
                else{
                    par1 -> par = par2;
                    if (par1 -> size == par2 -> size){
                        par2 -> size += 1;
                    }
                }
            }
        }
    }
    return;
}


void find_components(Node* graph, int w, int h) {
    int i, j;
    for (i= 0; i < w; i++) {
        for (j = 0; j < h; j++) {
            Node* node = &graph[i * w + j];
            union_set(node, node -> up);
            union_set(node, node -> down);
            union_set(node, node -> left);
            union_set(node, node -> right);
        }
    }
    return;
}


void color_components(Node* graph, int w, int h) {
    int i;
    unsigned char* output_image = malloc(w * h * 4 * sizeof(unsigned char));
    for (i = 0; i < w * h; i++) {
        Node* par = find_set(&graph[i]);
        if ((&graph[i]) -> R > 200){
            color_pixel(output_image, w, h, i, 250, 250, 250);
        }
        else{
            if (par == &graph[i]) {
                if (par -> size < 2) {
                    par -> R = 0;
                    par -> G = 0;
                    par -> B = 0;
                } else {
                    par -> R = rand() % 256;
                    par -> G = rand() % 256;
                    par -> B = rand() % 256;
                }
            }
            color_pixel(output_image, w, h, i, par -> R, par -> G, par -> B);
        }
    }

    char* filename = "output.png";
    lodepng_encode32_file(filename, output_image, w, h);
    free(output_image);
    return;
}


int main() {
    char *filename = "skull.png";
    int w, h, i, j;
    double Sigma = 0.95;
    double Gauss_matrix [5][5];
    make_Gauss_matrix(Gauss_matrix, Sigma);
    char *picture = load_png_file(filename, &w, &h);
    printf("w = %d   h = %d\n", w, h);
    Node* graph = malloc(w * h * sizeof(Node));
    int y = 0, x = 0;

    make_picture_grey(picture, w, h);

    Gauss_filter(picture, w, h, Sigma);

    applySobelFilter(picture, w, h);

    Gauss_filter(picture, w, h, 0.8);

    for (x = 0; x < w; x++) {
        for (y = 0; y < h; y++) {

            Node* node = &graph[y * w + x];


            unsigned char* pixel = picture[(y * w + x) * 4];
            node -> R = pixel + 0;
            node -> G = pixel + 1;
            node -> B = pixel + 2;
            node -> alpha= pixel + 3;

            if (y > 0) node -> up = &graph[(y - 1) * w + x];
            else node->up = NULL;

            if (y < h - 1 ) node -> down = &graph[(y + 1) * w + x];
            else node -> down = NULL;

            if (x > 0) node -> left = &graph[y * w + (x - 1)];
            else node -> left = NULL;

            if (x < w - 1) node -> right = &graph[y * w + (x + 1)];
            else node -> right = NULL;

            node -> par = node;
            node -> size = 0;
        }
    }
    find_components(graph, w, h);
    color_components(graph, w, h);
    free(graph);
    free(picture);

    return 0;
}
