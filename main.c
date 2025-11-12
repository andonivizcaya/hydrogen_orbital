#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "hydrogen.h"
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <stdio.h>

#define CLAY_IMPLEMENTATION
#include "clay.h"
#include "clay_renderer_raylib.c"

#define FPS 69
#define SCREEN_WIDTH 800.0f
#define SCREEN_HEIGHT 450.0f

#define MIN_RADIUS -50e-10
#define MAX_RADIUS 50e-10
#define MAX_VECTOR_SIZE 30

#define KNOB_RADIUS 16
#define MAX_STATES 8

#define BACKGROUND    CLITERAL (Color) { 18, 18, 18, 255 }
#define KNOB_NORMAL   CLITERAL (Color) { 255, 138, 50, 255 }
#define KNOB_DRAGGING CLITERAL (Color) { 225, 138, 50, 255 }

const Clay_Color COLOR_LIGHT          = (Clay_Color) { 224, 215, 210, 255 };
const Clay_Color COLOR_RED            = (Clay_Color) { 168, 66, 28, 255 };
const Clay_Color COLOR_ORANGE         = (Clay_Color) { 225, 138, 50, 255 };
const Clay_Color COLOR_BACKGROUND     = (Clay_Color) { 18, 18, 18, 255 };
const Clay_Color COLOR_TANSPARENT     = (Clay_Color) { 0, 0, 0, 0 };
const Clay_Color COLOR_KNOB_CONTAINER = (Clay_Color) { 50, 50, 50, 255 };
const Clay_Color COLOR_INFO           = (Clay_Color) { 30, 30, 30, 255 };

const uint32_t FONT_ID_BODY_24 = 0;

typedef struct {
    int *value;
    int min_value;
    int max_value;
    bool is_dragging;
} SliderState;

typedef enum {
    POINTS_MODE,
    MESH_MODE,
    POINTS_AND_MESH_MODE,
    WIREFRAME_MODE,
} RenderingMode;

size_t factorial(int n)
{
    if (n == 0) return 1;
    uint64_t res = 1;
    for (int i = 1; i <= n; ++i) res *= i;
    return res;
}

size_t doublefactorial(size_t n)
{
    if (n == 0 || n == 1) return 1;
    return n*doublefactorial(n - 2);
}

void hydrogen_matrix_compute_min_and_max_values(HydrogenMatrix *a, double min_value, double max_value)
{
    for (size_t i = 0; i < a->count; i++) {
        for (size_t j = 0; j < a->items[0].count; j++) {
            if (a->items[i].items[j] > max_value) max_value = a->items[i].items[j];
            if (a->items[i].items[j] < min_value) min_value = a->items[i].items[j];
        }
    }
}

void hydrogen_matrix_multiplication(HydrogenMatrix *a, HydrogenMatrix *b, HydrogenMatrix *c)
{
    if (a->items[0].count != b->count) {
        fprintf(stderr, "ERROR: columns of matrix `a` should be the samne as rows of matrix `b`. You passed %zu and %zu respectively\n", a[0].count, b->count);
        exit(1);
    }

    HydrogenRow c_row = {0};

    for (size_t i = 0; i < a->count; i++) {
        HydrogenRow c_row = {0};
        for (size_t j = 0; j < b->items[0].count; j++) {
            double result = 0.0L;
            for (size_t k = 0; k < b->count; k++) {
                result += a->items[i].items[k]*b->items[k].items[j];
            }
            da_append(&c_row, result);
        }
        da_append(c, c_row);
    }

    da_free(c_row);
}

void hydrogen_matrix_likewise_multiplication(HydrogenMatrix *a, HydrogenMatrix *b, HydrogenMatrix *c)
{
    if (a->items[0].count != b->items[0].count || a->count != b->count) {
        fprintf(stderr, "ERROR: matrices must have the same shape. [Matrix A]: %zux%zu, [Matrix B]: %zux%zu\n", a->count, a->items[0].count, b->count, b->items[0].count);
        exit(1);
    }

    HydrogenRow c_row = {0};

    for (size_t i = 0; i < a->count; i++) {
        HydrogenRow c_row = {0};
        for (size_t j = 0; j < b->items[0].count; j++) {
            da_append(&c_row, a->items[i].items[j]*b->items[i].items[j]);
        }
        da_append(c, c_row);
    }

    da_free(c_row);
}

double laguerre_polynomials(int n, double x)
{
    n = n - 1;
    if (fabs(x) > 1.0L) return 0.0L;
    if (n <= 0)         return 1.0L;
    if (n == 1)         return -x + 1.0L;

    return (((double)2.0L*n + 1.0L - x)*laguerre_polynomials(n, x) - (double)n*laguerre_polynomials(n - 1, x))/((double)(n + 1.0L));
}

double associated_laguerre_polynomials(int n, int l, double x)
{
    if (l == 0) return laguerre_polynomials(n, x);
    if (n <= 0) return 1.0L;
    if (n == 1) return 1.0L + l - x;
    return ((2.0L*n - 1.0L + l -x)*associated_laguerre_polynomials(n - 1, l, x) - (n - 1 + l)*associated_laguerre_polynomials(n - 2, l, x))/n;
}

double associated_legendre_function(int m, int l, double x)
{
    if (fabs(x) > 1.0L) return 0.0L;
    if (abs(m) > l)     return 0.0L;

    if (m == 0) {
        if (l == 0) return 1.0L;
        if (l == 1) return x;
        return (x*(double)(2.0L*l - 1)*associated_legendre_function(m, l - 1, x) - (double)(l + m - 1)*associated_legendre_function(m, l - 2, x))/(double)(l + m);
    }

    double mfactor = 1.0L;

    if (m < 0) {
        m = -1*m;
        mfactor = powf(-1.0L, (double)m)*((double)factorial(l - m)/(double)factorial(l + m));
    }

    if (m == l) {
        if (m == 0) return 1.0L;
        return mfactor*powf(-1.0, (double)l)*(double)doublefactorial(2*l - 1)*powf(1.0L - x*x, (double)l/2.0L);
    }

    if (m == l - 1) return x*(2.0L*(double)m + 1.0L)*associated_legendre_function(m, m, x);

    return mfactor*(x*(double)(2*l - 1)*associated_legendre_function(m, l - 1, x) - (double)(l + m - 1)*associated_legendre_function(m, l - 2, x))/(double)(l + m);
}

void hydrogen_matrix_format(HydrogenMatrix *matrix, bool with_size)
{
    printf("[\n");
    da_foreach(HydrogenRow, row, matrix) {
        printf("    [");
        da_foreach(double, value, row) {
            printf("%.10lf, ", *value);
        }
        printf("],");
        if (with_size) printf(" size: %zu", row->count);
        printf("\n");
    }
    printf("]\n");
    if (with_size) printf("%zux%zu\n", matrix->count, matrix->items[0].count);
}

void hydrogen_matrix_generate_wave_equation(HydrogenMatrix *spherical_normals, HydrogenMatrix *spherical_harmonics ,HydrogenMatrix *xs, HydrogenMatrix *ys, HydrogenMatrix *zs, size_t vector_size, size_t n, size_t l, int m)
{
    HydrogenMatrix cos_theta  = {0};
    HydrogenRow cos_theta_row = {0};
    HydrogenMatrix cos_phi    = {0};
    HydrogenRow cos_phi_row   = {0};

    HydrogenMatrix sin_theta  = {0};
    HydrogenRow sin_theta_row = {0};
    HydrogenMatrix sin_phi    = {0};
    HydrogenRow sin_phi_row   = {0};

    HydrogenMatrix ones  = {0};
    HydrogenRow ones_row = {0};

    HydrogenMatrix m_cos_phi  = {0};
    HydrogenRow m_cos_phi_row = {0};

    HydrogenMatrix legendre_polinomyals  = {0};
    HydrogenRow legendre_polinomyals_row = {0};

    HydrogenMatrix x_proy        = {0};
    HydrogenMatrix y_proy        = {0};
    HydrogenMatrix z_proy        = {0};
    HydrogenMatrix wave_equation = {0};

    for (size_t i = 0; i < 2*vector_size + 1; i++) {
        HydrogenRow cos_theta_row = {0};
        HydrogenRow sin_theta_row = {0};
        da_append(&cos_theta_row, (double)cos((double)i*M_PI/((double)2*vector_size)));
        da_append(&cos_theta, cos_theta_row);
        da_append(&sin_theta_row, (double)sin((double)i*M_PI/((double)2*vector_size)));
        da_append(&sin_theta, sin_theta_row);
    }

    for (int i = -1*(int)vector_size; i < (int)vector_size + 1; i++) {
        da_append(&cos_phi_row, (double)cos((double)i*M_PI/(double)vector_size));
        da_append(&sin_phi_row, (double)sin((double)i*M_PI/(double)vector_size));
        da_append(&ones_row, 1.0L);
        da_append(&m_cos_phi_row, (double)cos(abs(m)*(double)i*M_PI/(double)vector_size));
    }

    da_append(&cos_phi, cos_phi_row);
    da_append(&sin_phi, sin_phi_row);
    da_append(&ones, ones_row);
    da_append(&m_cos_phi, m_cos_phi_row);

    da_foreach(HydrogenRow, row, &cos_theta) {
        HydrogenRow legendre_polinomyals_row = {0};
        da_foreach(double, value, row) {
            da_append(&legendre_polinomyals_row, associated_legendre_function(m, l, *value));
        }
        da_append(&legendre_polinomyals, legendre_polinomyals_row);
    }

    hydrogen_matrix_multiplication(&legendre_polinomyals, &m_cos_phi, spherical_normals);
    hydrogen_matrix_multiplication(&legendre_polinomyals, &m_cos_phi, spherical_harmonics);
    da_foreach(HydrogenRow, row, spherical_harmonics) {
        da_foreach(double, value, row) *value = powf(fabs(*value), 1.0L/(l + 1.0L));
    }

    double step   = (MAX_RADIUS - MIN_RADIUS)/(spherical_harmonics->count*spherical_harmonics->items[0].count);
    double curr_x = MIN_RADIUS - step;

    HydrogenMatrix radial_part  = {0};
    HydrogenRow radial_part_row = {0};
    for (size_t i = 0; i < spherical_harmonics->count; i++) {
        HydrogenRow radial_part_row = {0};
        for (size_t j = 0; j < spherical_harmonics->items[0].count; j++) {
            curr_x += step;
            da_append(&radial_part_row, associated_laguerre_polynomials(n - l - 1, 2*l + 1, curr_x));
        }
        da_append(&radial_part, radial_part_row);
    }

    hydrogen_matrix_multiplication(&sin_theta, &cos_phi, &x_proy);
    hydrogen_matrix_multiplication(&sin_theta, &sin_phi, &y_proy);
    hydrogen_matrix_multiplication(&cos_theta, &ones, &z_proy);

    hydrogen_matrix_likewise_multiplication(&radial_part, spherical_harmonics, &wave_equation);
    hydrogen_matrix_likewise_multiplication(&wave_equation, &x_proy, xs);
    hydrogen_matrix_likewise_multiplication(&wave_equation, &y_proy, ys);
    hydrogen_matrix_likewise_multiplication(&wave_equation, &z_proy, zs);

    da_free(cos_theta_row);
    da_free(sin_theta_row);
    da_free(cos_phi_row);
    da_free(sin_phi_row);
    da_free(m_cos_phi_row);
    da_free(m_cos_phi);
    da_free(sin_theta);
    da_free(cos_theta);
    da_free(sin_phi);
    da_free(cos_phi);
    da_free(ones_row);
    da_free(ones);
    da_free(legendre_polinomyals_row);
    da_free(legendre_polinomyals);
    da_free(radial_part_row);
    da_free(radial_part);
    da_free(x_proy);
    da_free(y_proy);
    da_free(z_proy);
    da_free(wave_equation);
}

SliderState slider_n = {0};
SliderState slider_l = {0};
SliderState slider_m = {0};

void handle_clay_errors(Clay_ErrorData errorData)
{
    printf("%s", errorData.errorText.chars);
    switch(errorData.errorType) {
        case CLAY_ERROR_TYPE_TEXT_MEASUREMENT_FUNCTION_NOT_PROVIDED: break;
        case CLAY_ERROR_TYPE_ARENA_CAPACITY_EXCEEDED: break;
        case CLAY_ERROR_TYPE_ELEMENTS_CAPACITY_EXCEEDED: break;
        case CLAY_ERROR_TYPE_TEXT_MEASUREMENT_CAPACITY_EXCEEDED: break;
        case CLAY_ERROR_TYPE_DUPLICATE_ID: break;
        case CLAY_ERROR_TYPE_FLOATING_CONTAINER_PARENT_NOT_FOUND: break;
        case CLAY_ERROR_TYPE_PERCENTAGE_OVER_1: break;
        case CLAY_ERROR_TYPE_INTERNAL_ERROR: break;

        default: break;
    }
}

Color get_color(float v,float vmin, float vmax)
{
    Color c = { 255, 255, 255, 255 };
    float dv;

    if (v < vmin) v = vmin;
    if (v > vmax) v = vmax;
    dv = vmax - vmin;

    if (v < (vmin + 0.25*dv)) {
       c.r = 0;
       c.g = (unsigned char)255*(4*(v - vmin)/dv);
    } else if (v < (vmin + 0.5*dv)) {
       c.r = 0;
       c.b = (unsigned char)255*(1 + 4*(vmin + 0.25*dv - v)/dv);
    } else if (v < (vmin + 0.75*dv)) {
       c.r = (unsigned char)255*4*(v - vmin - 0.5*dv)/dv;
       c.b = 0;
    } else {
       c.g = (unsigned char)255*(1 + 4*(vmin + 0.75*dv - v)/dv);
       c.b = 0;
    }

    return c;
}

Mesh generate_mesh_from_points(HydrogenMatrix *xs, HydrogenMatrix *ys, HydrogenMatrix *zs, HydrogenMatrix *values)
{
    size_t rows = xs->count;
    size_t cols = xs->items[0].count;

    size_t vertex_count   = rows*cols;
    size_t triangle_count = (rows - 1)*(cols - 1)*2;

    Mesh mesh          = { 0 };
    mesh.vertexCount   = vertex_count;
    mesh.triangleCount = triangle_count;

    mesh.vertices  = (float *)MemAlloc(vertex_count*3*sizeof(float));
    mesh.normals   = (float *)MemAlloc(vertex_count*3*sizeof(float));
    mesh.texcoords = (float *)MemAlloc(vertex_count*2*sizeof(float));
    mesh.colors    = (unsigned char *)MemAlloc(vertex_count*4*sizeof(unsigned char));
    mesh.indices   = (unsigned short *)MemAlloc(triangle_count*3*sizeof(unsigned short));

    float min_val = 1e10, max_val = -1e10;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            float val = fabs(values->items[i].items[j]);
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            size_t idx = i*cols + j;

            mesh.vertices[idx*3 + 0] = (float)xs->items[i].items[j];
            mesh.vertices[idx*3 + 1] = (float)ys->items[i].items[j];
            mesh.vertices[idx*3 + 2] = (float)zs->items[i].items[j];

            float magnitude = fabs(values->items[i].items[j]);
            Color jet_color = get_color(magnitude, min_val, max_val);
            mesh.colors[idx*4 + 0] = jet_color.r;
            mesh.colors[idx*4 + 1] = jet_color.g;
            mesh.colors[idx*4 + 2] = jet_color.b;
            mesh.colors[idx*4 + 3] = 200;

            mesh.texcoords[idx*2 + 0] = (float)j/(float)cols;
            mesh.texcoords[idx*2 + 1] = (float)i/(float)rows;

            mesh.normals[idx*3 + 0] = 0.0f;
            mesh.normals[idx*3 + 1] = 0.0f;
            mesh.normals[idx*3 + 2] = 0.0f;
        }
    }

    size_t tri_idx = 0;
    for (size_t i = 0; i < rows - 1; i++) {
        for (size_t j = 0; j < cols - 1; j++) {
            size_t v0 = i*cols + j;
            size_t v1 = i*cols + (j + 1);
            size_t v2 = (i + 1)*cols + j;
            size_t v3 = (i + 1)*cols + (j + 1);

            mesh.indices[tri_idx*3 + 0] = v0;
            mesh.indices[tri_idx*3 + 1] = v1;
            mesh.indices[tri_idx*3 + 2] = v2;
            tri_idx++;

            mesh.indices[tri_idx*3 + 0] = v1;
            mesh.indices[tri_idx*3 + 1] = v3;
            mesh.indices[tri_idx*3 + 2] = v2;
            tri_idx++;
        }
    }

    mesh.triangleCount = tri_idx;

    for (int i = 0; i < mesh.triangleCount; i++) {
        unsigned short idx0 = mesh.indices[i*3 + 0];
        unsigned short idx1 = mesh.indices[i*3 + 1];
        unsigned short idx2 = mesh.indices[i*3 + 2];

        Vector3 v0 = { mesh.vertices[idx0*3], mesh.vertices[idx0*3 + 1], mesh.vertices[idx0*3 + 2] };
        Vector3 v1 = { mesh.vertices[idx1*3], mesh.vertices[idx1*3 + 1], mesh.vertices[idx1*3 + 2] };
        Vector3 v2 = { mesh.vertices[idx2*3], mesh.vertices[idx2*3 + 1], mesh.vertices[idx2*3 + 2] };

        Vector3 edge1  = Vector3Subtract(v1, v0);
        Vector3 edge2  = Vector3Subtract(v2, v0);
        Vector3 normal = Vector3CrossProduct(edge1, edge2);

        mesh.normals[idx0*3 + 0] += normal.x;
        mesh.normals[idx0*3 + 1] += normal.y;
        mesh.normals[idx0*3 + 2] += normal.z;

        mesh.normals[idx1*3 + 0] += normal.x;
        mesh.normals[idx1*3 + 1] += normal.y;
        mesh.normals[idx1*3 + 2] += normal.z;

        mesh.normals[idx2*3 + 0] += normal.x;
        mesh.normals[idx2*3 + 1] += normal.y;
        mesh.normals[idx2*3 + 2] += normal.z;
    }

    for (size_t i = 0; i < vertex_count; i++) {
        Vector3 normal = {
            mesh.normals[i*3 + 0],
            mesh.normals[i*3 + 1],
            mesh.normals[i*3 + 2]
        };

        float len = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
        if (len > 0.0001f) {
            normal = Vector3Normalize(normal);
            mesh.normals[i*3 + 0] = normal.x;
            mesh.normals[i*3 + 1] = normal.y;
            mesh.normals[i*3 + 2] = normal.z;
        } else {
            Vector3 pos = {mesh.vertices[i*3], mesh.vertices[i*3 + 1], mesh.vertices[i*3 + 2]};
            float r = sqrtf(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
            if (r > 0.0001f) {
                mesh.normals[i*3 + 0] = pos.x / r;
                mesh.normals[i*3 + 1] = pos.y / r;
                mesh.normals[i*3 + 2] = pos.z / r;
            } else {
                mesh.normals[i*3 + 0] = 0.0f;
                mesh.normals[i*3 + 1] = 1.0f;
                mesh.normals[i*3 + 2] = 0.0f;
            }
        }
    }

    UploadMesh(&mesh, false);

    return mesh;
}

void slider_handle_interaction(Clay_String slider_id, SliderState *slider_state, Vector2 mouse_position, bool is_mouse_down)
{
    Clay_ElementId element_id = Clay_GetElementId(slider_id);
    Clay_ElementData element  = Clay_GetElementData(element_id);

    if (!element.found) return;

    bool is_hovering = Clay_PointerOver(element_id);

    if (is_hovering && is_mouse_down) slider_state->is_dragging = true;
    if (!is_mouse_down)               slider_state->is_dragging = false;

    if (slider_state->is_dragging) {
        float bar_x      = element.boundingBox.x;
        float bar_width  = element.boundingBox.width;
        float relative_x = mouse_position.x - bar_x;

        relative_x = fmax(0.0f, fmin(relative_x, bar_width));

        float normalized = relative_x/bar_width;
        int new_value    = slider_state->min_value + (int)(normalized*(slider_state->max_value - slider_state->min_value));

        *slider_state->value = new_value;
    }
}

float slider_get_knob_position(SliderState *slider_state, float bar_width)
{
    if (slider_state->max_value == slider_state->min_value) return 0.0f;

    float normalized = (float)(*slider_state->value - slider_state->min_value)/(float)(slider_state->max_value - slider_state->min_value);
    return normalized*bar_width;
}

int main(void)
{
    Clay_Raylib_Initialize((int)SCREEN_WIDTH, (int)SCREEN_HEIGHT, "Hydrogen Atom GIGACHAD", FLAG_BORDERLESS_WINDOWED_MODE | FLAG_VSYNC_HINT | FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);

    uint64_t totalMemorySize = Clay_MinMemorySize();
    Clay_Arena arena         = Clay_CreateArenaWithCapacityAndMemory(totalMemorySize, malloc(totalMemorySize));

    Clay_Initialize(arena, (Clay_Dimensions) { SCREEN_WIDTH, SCREEN_HEIGHT }, (Clay_ErrorHandler) { handle_clay_errors, 0 });
    Font fonts[2];
    fonts[FONT_ID_BODY_24] = LoadFontEx("resources/fonts/Iosevka-Regular.ttc", 48, 0, 400);

    Clay_SetMeasureTextFunction(Raylib_MeasureText, fonts);

    SetTargetFPS(FPS);

    Camera3D camera   = {0};
    camera.position   = (Vector3){ 5.0f, 5.0f, 5.0f };
    camera.target     = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up         = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy       = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    int n = 1;
    int l = 0;
    int m = 0;

    slider_n.value       = &n;
    slider_n.min_value   = 1;
    slider_n.max_value   = 7;
    slider_n.is_dragging = false;

    slider_l.value       = &l;
    slider_l.min_value   = 0;
    slider_l.max_value   = n - 1;
    slider_l.is_dragging = false;

    slider_m.value       = &m;
    slider_m.min_value   = -l;
    slider_m.max_value   = l;
    slider_m.is_dragging = false;

    size_t vector_size = max(MAX_VECTOR_SIZE, 6*l);

    HydrogenMatrix spherical_harmonics = {0};
    HydrogenMatrix spherical_normals   = {0};
    HydrogenMatrix xs                  = {0};
    HydrogenMatrix ys                  = {0};
    HydrogenMatrix zs                  = {0};

    hydrogen_matrix_generate_wave_equation(&spherical_normals, &spherical_harmonics, &xs, &ys, &zs, vector_size, (size_t)n, (size_t)l, m);

    Mesh orbital_mesh = generate_mesh_from_points(&xs, &ys, &zs, &spherical_normals);
    Material material = LoadMaterialDefault();

    material.maps[MATERIAL_MAP_ALBEDO].color = WHITE;
    material.maps[MATERIAL_MAP_EMISSION].color = WHITE;
    material.maps[MATERIAL_MAP_EMISSION].value = 1.0f;

    bool needs_regeneration      = false;
    bool has_mesh                = false;
    RenderingMode rendering_mode = POINTS_MODE;
    bool menu_collapsed          = false;
    bool is_cursor_on_menu       = false;
    bool show_axis               = true;

    while (!WindowShouldClose()) {
        int old_n = n;
        int old_l = l;
        int old_m = m;

        Clay_SetLayoutDimensions((Clay_Dimensions) { SCREEN_WIDTH, SCREEN_HEIGHT });

        Vector2 mouse_position   = GetMousePosition();
        Vector2 scroll_delta     = GetMouseWheelMoveV();
        bool is_left_mouse_down  = IsMouseButtonDown(MOUSE_LEFT_BUTTON);
        bool is_right_mouse_down = IsMouseButtonDown(MOUSE_RIGHT_BUTTON);

        Clay_SetPointerState((Clay_Vector2) { mouse_position.x, mouse_position.y }, is_left_mouse_down);
        Clay_UpdateScrollContainers(true, (Clay_Vector2) { scroll_delta.x, scroll_delta.y }, GetFrameTime());

        if (!is_right_mouse_down && !menu_collapsed) {
            slider_handle_interaction(CLAY_STRING("SliderBarN"), &slider_n, mouse_position, is_left_mouse_down);
            slider_handle_interaction(CLAY_STRING("SliderBarL"), &slider_l, mouse_position, is_left_mouse_down);
            slider_handle_interaction(CLAY_STRING("SliderBarM"), &slider_m, mouse_position, is_left_mouse_down);
        }

        if (l >= n) l = n - 1;
        if (l < 0)  l = 0;
        slider_l.max_value = n - 1;

        if (m > l)  m = l;
        if (m < -l) m = -l;
        slider_m.min_value = -l;
        slider_m.max_value = l;

        if (old_n != n || old_l != l || old_m != m) needs_regeneration = true;

        if (needs_regeneration) {
            if (has_mesh) UnloadMesh(orbital_mesh);

            da_free(spherical_harmonics);
            da_free(spherical_normals);
            da_free(xs);
            da_free(ys);
            da_free(zs);

            spherical_harmonics = (HydrogenMatrix){0};
            spherical_normals   = (HydrogenMatrix){0};

            xs = (HydrogenMatrix){0};
            ys = (HydrogenMatrix){0};
            zs = (HydrogenMatrix){0};

            vector_size = max(MAX_VECTOR_SIZE, 6*l);
            hydrogen_matrix_generate_wave_equation(&spherical_normals, &spherical_harmonics, &xs, &ys, &zs, vector_size, (size_t)n, (size_t)l, m);

            orbital_mesh = generate_mesh_from_points(&xs, &ys, &zs, &spherical_normals);
            has_mesh     = true;

            needs_regeneration = false;
        }

        BeginDrawing();
        ClearBackground(BACKGROUND);

        if (!slider_n.is_dragging && !slider_l.is_dragging && !slider_m.is_dragging) {
            UpdateCameraPro(&camera,
                (Vector3){
                    (IsKeyDown(KEY_W) || IsKeyDown(KEY_UP) || (is_right_mouse_down && is_left_mouse_down))*0.1f - (IsKeyDown(KEY_S) || IsKeyDown(KEY_DOWN))*0.1f,
                    (IsKeyDown(KEY_D) || IsKeyDown(KEY_RIGHT))*0.1f                                             - (IsKeyDown(KEY_A) || IsKeyDown(KEY_LEFT))*0.1f,
                    IsKeyDown(KEY_SPACE)*0.1f - IsKeyDown(KEY_X)*0.1f
                },
                (Vector3){
                    (is_left_mouse_down || is_right_mouse_down) && !is_cursor_on_menu ? -GetMouseDelta().x*0.1f : 0.0f,
                    (is_left_mouse_down || is_right_mouse_down) && !is_cursor_on_menu ? -GetMouseDelta().y*0.1f : 0.0f,
                    0.0f
                },
                !is_cursor_on_menu ? GetMouseWheelMove()*2.0f : 0.0f
            );
        }

        if (IsKeyPressed(KEY_ONE)) {
            has_mesh = false;
            rendering_mode = POINTS_MODE;
        }
        if (IsKeyPressed(KEY_TWO)) {
            has_mesh = true;
            rendering_mode = MESH_MODE;
        }
        if (IsKeyPressed(KEY_THREE)) {
            has_mesh = true;
            rendering_mode = POINTS_AND_MESH_MODE;
        }
        if (IsKeyPressed(KEY_FOUR)) {
            has_mesh = true;
            rendering_mode = WIREFRAME_MODE;
        }
        if (IsKeyPressed(KEY_TAB)) {
            menu_collapsed = !menu_collapsed;
        }

        if (IsKeyPressed(KEY_O)) {
            show_axis = !show_axis;
        }

        if (IsKeyPressed(KEY_R)) {
            camera.position   = (Vector3){ 5.0f, 5.0f, 5.0f };
            camera.target     = (Vector3){ 0.0f, 0.0f, 0.0f };
            camera.up         = (Vector3){ 0.0f, 1.0f, 0.0f };
            camera.fovy       = 45.0f;
        }

        BeginMode3D(camera);
        switch (rendering_mode) {
            case POINTS_MODE: {
                for (size_t i = 0; i < spherical_harmonics.count; i++) {
                    for (size_t j = 0; j < spherical_harmonics.items[0].count; j++) {
                        DrawPoint3D((Vector3){ (float)xs.items[i].items[j], (float)ys.items[i].items[j], (float)zs.items[i].items[j] }, RED);
                    }
                }
                break;
            }
            case MESH_MODE: {
                if (has_mesh) {
                    rlDisableBackfaceCulling();
                    rlBegin(RL_TRIANGLES);
                    for (int i = 0; i < orbital_mesh.triangleCount; i++) {
                        unsigned short idx0 = orbital_mesh.indices[i*3 + 0];
                        unsigned short idx1 = orbital_mesh.indices[i*3 + 1];
                        unsigned short idx2 = orbital_mesh.indices[i*3 + 2];

                        rlColor4ub(orbital_mesh.colors[idx0*4], orbital_mesh.colors[idx0*4 + 1], orbital_mesh.colors[idx0*4 + 2], orbital_mesh.colors[idx0*4 + 3]);
                        rlVertex3f(orbital_mesh.vertices[idx0*3], orbital_mesh.vertices[idx0*3 + 1], orbital_mesh.vertices[idx0*3 + 2]);

                        rlColor4ub(orbital_mesh.colors[idx1*4], orbital_mesh.colors[idx1*4 + 1], orbital_mesh.colors[idx1*4 + 2], orbital_mesh.colors[idx1*4 + 3]);
                        rlVertex3f(orbital_mesh.vertices[idx1*3], orbital_mesh.vertices[idx1*3 + 1], orbital_mesh.vertices[idx1*3 + 2]);

                        rlColor4ub(orbital_mesh.colors[idx2*4], orbital_mesh.colors[idx2*4 + 1], orbital_mesh.colors[idx2*4 + 2], orbital_mesh.colors[idx2*4 + 3]);
                        rlVertex3f(orbital_mesh.vertices[idx2*3], orbital_mesh.vertices[idx2*3 + 1], orbital_mesh.vertices[idx2*3 + 2]);
                    }
                    rlEnd();
                    rlEnableBackfaceCulling();
                }
                break;
            }
            case POINTS_AND_MESH_MODE: {
                if (has_mesh) {
                    rlDisableBackfaceCulling();
                    rlBegin(RL_TRIANGLES);
                    for (int i = 0; i < orbital_mesh.triangleCount; i++) {
                        unsigned short idx0 = orbital_mesh.indices[i*3 + 0];
                        unsigned short idx1 = orbital_mesh.indices[i*3 + 1];
                        unsigned short idx2 = orbital_mesh.indices[i*3 + 2];

                        rlColor4ub(orbital_mesh.colors[idx0*4], orbital_mesh.colors[idx0*4 + 1], orbital_mesh.colors[idx0*4 + 2], orbital_mesh.colors[idx0*4 + 3]);
                        rlVertex3f(orbital_mesh.vertices[idx0*3], orbital_mesh.vertices[idx0*3 + 1], orbital_mesh.vertices[idx0*3 + 2]);

                        rlColor4ub(orbital_mesh.colors[idx1*4], orbital_mesh.colors[idx1*4 + 1], orbital_mesh.colors[idx1*4 + 2], orbital_mesh.colors[idx1*4 + 3]);
                        rlVertex3f(orbital_mesh.vertices[idx1*3], orbital_mesh.vertices[idx1*3 + 1], orbital_mesh.vertices[idx1*3 + 2]);

                        rlColor4ub(orbital_mesh.colors[idx2*4], orbital_mesh.colors[idx2*4 + 1], orbital_mesh.colors[idx2*4 + 2], orbital_mesh.colors[idx2*4 + 3]);
                        rlVertex3f(orbital_mesh.vertices[idx2*3], orbital_mesh.vertices[idx2*3 + 1], orbital_mesh.vertices[idx2*3 + 2]);
                    }
                    rlEnd();
                    rlEnableBackfaceCulling();
                }
                for (size_t i = 0; i < spherical_harmonics.count; i++) {
                    for (size_t j = 0; j < spherical_harmonics.items[0].count; j++) {
                        DrawPoint3D((Vector3){ (float)xs.items[i].items[j], (float)ys.items[i].items[j], (float)zs.items[i].items[j] }, RED);
                    }
                }
                break;
            }
            case WIREFRAME_MODE: {
                if (has_mesh) {
                    for (int i = 0; i < orbital_mesh.triangleCount; i++) {
                        unsigned short idx0 = orbital_mesh.indices[i*3 + 0];
                        unsigned short idx1 = orbital_mesh.indices[i*3 + 1];
                        unsigned short idx2 = orbital_mesh.indices[i*3 + 2];

                        Vector3 v0 = { orbital_mesh.vertices[idx0*3], orbital_mesh.vertices[idx0*3 + 1], orbital_mesh.vertices[idx0*3 + 2] };
                        Vector3 v1 = { orbital_mesh.vertices[idx1*3], orbital_mesh.vertices[idx1*3 + 1], orbital_mesh.vertices[idx1*3 + 2] };
                        Vector3 v2 = { orbital_mesh.vertices[idx2*3], orbital_mesh.vertices[idx2*3 + 1], orbital_mesh.vertices[idx2*3 + 2] };

                        DrawLine3D(v0, v1, ORANGE);
                        DrawLine3D(v1, v2, ORANGE);
                        DrawLine3D(v2, v0, ORANGE);
                    }
                }
                break;
            }
            default: break;
        }

        if (show_axis) {
            DrawLine3D((Vector3){0, 0, 0}, (Vector3){1000, 0, 0}, RED);
            DrawLine3D((Vector3){0, 0, 0}, (Vector3){0, 1000, 0}, GREEN);
            DrawLine3D((Vector3){0, 0, 0}, (Vector3){0, 0, 1000}, BLUE);
        }

        EndMode3D();

        Clay_BeginLayout();

        CLAY(CLAY_ID("OuterContainer"), { .layout = { .sizing = { CLAY_SIZING_GROW(0), CLAY_SIZING_GROW(0) }, .padding = CLAY_PADDING_ALL(16), .childGap = 16 }, .backgroundColor = COLOR_TANSPARENT }) {
            if (menu_collapsed) {
                CLAY(CLAY_ID("CollapsedMenu"), {
                    .layout = {
                        .layoutDirection = CLAY_TOP_TO_BOTTOM,
                        .sizing = { .width = CLAY_SIZING_FIXED(60), .height = CLAY_SIZING_FIXED(60) },
                        .padding = CLAY_PADDING_ALL(8),
                        .childAlignment = { .x = CLAY_ALIGN_X_CENTER, .y = CLAY_ALIGN_Y_CENTER }
                    },
                    .backgroundColor = COLOR_BACKGROUND,
                    .cornerRadius = CLAY_CORNER_RADIUS(8)
                }) {
                    CLAY_TEXT(CLAY_STRING("☰"), CLAY_TEXT_CONFIG({ .fontSize = 32, .textColor = COLOR_ORANGE }));
                }
            } else {
                CLAY(CLAY_ID("SideBar"), {
                    .layout = {
                        .layoutDirection = CLAY_TOP_TO_BOTTOM,
                        .sizing = { .width = CLAY_SIZING_FIXED(300), .height = CLAY_SIZING_GROW(0) },
                        .padding = CLAY_PADDING_ALL(16),
                        .childGap = 12
                    },
                    .backgroundColor = COLOR_BACKGROUND,
                    .cornerRadius = CLAY_CORNER_RADIUS(8)
                }) {
                    CLAY(CLAY_ID("TitleBar"), {
                        .layout = {
                            .sizing = { .width = CLAY_SIZING_GROW(0), .height = CLAY_SIZING_GROW(0) },
                            .padding = CLAY_PADDING_ALL(8),
                            .childAlignment = { .x = CLAY_ALIGN_X_CENTER, .y = CLAY_ALIGN_Y_CENTER }
                        },
                        .backgroundColor = COLOR_RED,
                        .cornerRadius = CLAY_CORNER_RADIUS(6)
                    }) {
                        CLAY_TEXT(CLAY_STRING("Quantum Numbers [TAB]"), CLAY_TEXT_CONFIG({ .fontSize = 20, .textColor = COLOR_LIGHT }));
                    }

                    CLAY(CLAY_ID("SliderContainerN"), { .layout = { .layoutDirection = CLAY_TOP_TO_BOTTOM, .sizing = { .width = CLAY_SIZING_GROW(0) }, .childGap = 8 } }) {
                        char n_label[32];
                        snprintf(n_label, sizeof(n_label), "n = %d", n);
                        CLAY_TEXT(((Clay_String) { .length = strlen(n_label), .chars = n_label }), CLAY_TEXT_CONFIG({ .fontSize = 18, .textColor = COLOR_LIGHT }));

                        CLAY(CLAY_ID("SliderBarN"), { .layout = { .sizing = { .width = CLAY_SIZING_GROW(0), .height = CLAY_SIZING_FIXED(40) }, .padding = CLAY_PADDING_ALL(4), .childAlignment = { .x = CLAY_ALIGN_X_LEFT, .y = CLAY_ALIGN_Y_CENTER } }, .backgroundColor = COLOR_KNOB_CONTAINER, .cornerRadius = CLAY_CORNER_RADIUS(20) }) {}
                    }

                    CLAY(CLAY_ID("SliderContainerL"), { .layout = { .layoutDirection = CLAY_TOP_TO_BOTTOM, .sizing = { .width = CLAY_SIZING_GROW(0) }, .childGap = 8 } }) {
                        char l_label[32];
                        snprintf(l_label, sizeof(l_label), "l = %d", l);
                        CLAY_TEXT(((Clay_String) { .length = strlen(l_label), .chars = l_label }), CLAY_TEXT_CONFIG({ .fontSize = 18, .textColor = COLOR_LIGHT }));

                        CLAY(CLAY_ID("SliderBarL"), { .layout = { .sizing = { .width = CLAY_SIZING_GROW(0), .height = CLAY_SIZING_FIXED(40) }, .padding = CLAY_PADDING_ALL(4), .childAlignment = { .x = CLAY_ALIGN_X_LEFT, .y = CLAY_ALIGN_Y_CENTER } }, .backgroundColor = COLOR_KNOB_CONTAINER, .cornerRadius = CLAY_CORNER_RADIUS(20) }) {}
                    }

                    CLAY(CLAY_ID("SliderContainerM"), { .layout = { .layoutDirection = CLAY_TOP_TO_BOTTOM, .sizing = { .width = CLAY_SIZING_GROW(0) }, .childGap = 8 } }) {
                        char m_label[32];
                        snprintf(m_label, sizeof(m_label), "m = %d", m);
                        CLAY_TEXT(((Clay_String) { .length = strlen(m_label), .chars = m_label }), CLAY_TEXT_CONFIG({ .fontSize = 18, .textColor = COLOR_LIGHT }));

                        CLAY(CLAY_ID("SliderBarM"), { .layout = { .sizing = { .width = CLAY_SIZING_GROW(0), .height = CLAY_SIZING_FIXED(40) }, .padding = CLAY_PADDING_ALL(4), .childAlignment = { .x = CLAY_ALIGN_X_LEFT, .y = CLAY_ALIGN_Y_CENTER } }, .backgroundColor = COLOR_KNOB_CONTAINER, .cornerRadius = CLAY_CORNER_RADIUS(20) }) {}
                    }

                    CLAY(CLAY_ID("InfoSection"), {
                        .layout = {
                            .layoutDirection = CLAY_TOP_TO_BOTTOM,
                            .sizing = { .width = CLAY_SIZING_GROW(0), .height = CLAY_SIZING_GROW(0) },
                            .padding = CLAY_PADDING_ALL(12),
                            .childGap = 3
                        },
                        .backgroundColor = COLOR_INFO,
                        .cornerRadius = CLAY_CORNER_RADIUS(6),
                        .clip = { .vertical = true, .childOffset = Clay_GetScrollOffset() }
                    }) {
                        CLAY_TEXT(CLAY_STRING("Controls"), CLAY_TEXT_CONFIG({ .fontSize = 16, .textColor = COLOR_ORANGE }));
                        CLAY_TEXT(CLAY_STRING(""), CLAY_TEXT_CONFIG({ .fontSize = 6, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("Camera:"), CLAY_TEXT_CONFIG({ .fontSize = 13, .textColor = COLOR_ORANGE }));
                        CLAY_TEXT(CLAY_STRING("  WASD - Move"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  Mouse - Rotate"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  Scroll - Zoom"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING(""), CLAY_TEXT_CONFIG({ .fontSize = 6, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("Rendering:"), CLAY_TEXT_CONFIG({ .fontSize = 13, .textColor = COLOR_ORANGE }));
                        CLAY_TEXT(CLAY_STRING("  1 - Points"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  2 - Mesh"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  3 - Both"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  4 - Wireframe"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING(""), CLAY_TEXT_CONFIG({ .fontSize = 6, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("UI:"), CLAY_TEXT_CONFIG({ .fontSize = 13, .textColor = COLOR_ORANGE }));
                        CLAY_TEXT(CLAY_STRING("  TAB - Toggle menu"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  O   - Show axis"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                        CLAY_TEXT(CLAY_STRING("  R   - Reset camera"), CLAY_TEXT_CONFIG({ .fontSize = 11, .textColor = COLOR_LIGHT }));
                    }
                }
            }
        }

        Clay_RenderCommandArray commands = Clay_EndLayout();
        Clay_Raylib_Render(commands, fonts);

        if (!menu_collapsed) {
            if (Clay_PointerOver(Clay_GetElementId(CLAY_STRING("OuterContainer")))) is_cursor_on_menu = false;
            if (Clay_PointerOver(Clay_GetElementId(CLAY_STRING("SideBar"))))        is_cursor_on_menu = true;

            Clay_ElementData slider_bar_n = Clay_GetElementData(Clay_GetElementId(CLAY_STRING("SliderBarN")));
            if (slider_bar_n.found) {
                float knob_pos   = slider_get_knob_position(&slider_n, slider_bar_n.boundingBox.width - MAX_STATES);
                float knob_x     = slider_bar_n.boundingBox.x + KNOB_RADIUS/2 + knob_pos;
                float knob_y     = slider_bar_n.boundingBox.y + slider_bar_n.boundingBox.height/2;
                Color knob_color = slider_n.is_dragging ? KNOB_NORMAL : KNOB_DRAGGING;
                if (slider_n.max_value == n) knob_x = knob_x - KNOB_RADIUS;
                DrawCircle((int)knob_x, (int)knob_y, KNOB_RADIUS, knob_color);
            }

            Clay_ElementData slider_bar_l = Clay_GetElementData(Clay_GetElementId(CLAY_STRING("SliderBarL")));
            if (slider_bar_l.found) {
                float knob_pos   = slider_get_knob_position(&slider_l, slider_bar_l.boundingBox.width - MAX_STATES);
                float knob_x     = slider_bar_l.boundingBox.x + KNOB_RADIUS/2 + knob_pos;
                float knob_y     = slider_bar_l.boundingBox.y + slider_bar_l.boundingBox.height/2;
                Color knob_color = slider_l.is_dragging ? KNOB_NORMAL : KNOB_DRAGGING;
                if (slider_l.max_value == l && slider_l.min_value != slider_l.max_value) knob_x = knob_x - KNOB_RADIUS;
                DrawCircle((int)knob_x, (int)knob_y, KNOB_RADIUS, knob_color);
            }

            Clay_ElementData slider_bar_m = Clay_GetElementData(Clay_GetElementId(CLAY_STRING("SliderBarM")));
            if (slider_bar_m.found) {
                float knob_pos   = slider_get_knob_position(&slider_m, slider_bar_m.boundingBox.width - MAX_STATES);
                float knob_x     = slider_bar_m.boundingBox.x + KNOB_RADIUS/2 + knob_pos;
                float knob_y     = slider_bar_m.boundingBox.y + slider_bar_m.boundingBox.height/2;
                Color knob_color = slider_m.is_dragging ? KNOB_NORMAL : KNOB_DRAGGING;
                if (slider_m.max_value == m && slider_m.min_value != slider_m.max_value) knob_x = knob_x - KNOB_RADIUS;
                DrawCircle((int)knob_x, (int)knob_y, KNOB_RADIUS, knob_color);
            }
        } else {
            is_cursor_on_menu = false;
        }

        EndDrawing();
    }

    UnloadMesh(orbital_mesh);
    UnloadMaterial(material);
    CloseWindow();

    return 0;
}
