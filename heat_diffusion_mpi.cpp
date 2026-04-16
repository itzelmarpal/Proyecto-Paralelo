#include<iostream>
#include<vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include<mpi.h>
using namespace std;

int N = 50;
double ALPHA = 0.1;
int ITERATIONS = 2000;
int recs = 10;

double heat_diffusion(double u_ij,
                      double u_im1j, 
                      double u_ip1j,
                      double u_ijm1,
                      double u_ijp1) {
    return u_ij + ALPHA * (u_im1j + u_ip1j + u_ijm1 + u_ijp1 
        - 4 * u_ij);
}

void write_vtx(int nrec, const vector<double>& U) {
    // Filename
    stringstream ss;
    ss << "S-" << setw(4) << setfill('0') << nrec << ".vtk";
    string filename = ss.str();

    ofstream file(filename);
    file << "# vtk DataFile Version 3.0\n";
    file << "Heat diffusion\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    file << "\nDIMENSIONS " << N << " " << N << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING 1 1 1\n";

    file << "\nPOINT_DATA " << N * N << "\n";
    file << "SCALARS temperature float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            file << scientific << setprecision(5) << U[i * N + j] << "\n";
        }
    }

    file.close();

    cout << "File " << nrec << " created" << endl;
}


int main() {
    int n_tasks, task_id;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    // ----- Domain divison -----
    // Width of processes 0, 1, ..., (n_tasks - 2)
    int width = (N - 2) / n_tasks;  
    // Width of process (n_tasks - 1)
    int width_f = (N - 1) - width * (n_tasks - 1);  
    // Width of current process
    int w = (task_id == n_tasks - 1 ? width_f : width);
    int start_j = task_id * width + 1;
    vector<double> u(N * (w + 2), 0);
    vector<double> u_new(N * (w + 2), 0);

    // Check if  (N/2, N/2) is contained
    bool is_center_contained = (start_j <= N/2) && (N/2 < start_j + w);
    if (is_center_contained) {
        int local_j = N / 2 - start_j + 1;
        u[N/2 * (w+2) + local_j] = 100.0;
    }

    // Local column type (for u)
    MPI_Datatype column;
    MPI_Type_vector(N, 1, w + 2, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    // Global block of N x width (for u_global)
    MPI_Datatype global_block;
    MPI_Type_vector(N, width, N, MPI_DOUBLE, &global_block);
    MPI_Type_commit(&global_block);

    // Global block of N x width_f (for u_global)
    MPI_Datatype global_block_f;
    MPI_Type_vector(N, width_f, N, MPI_DOUBLE, &global_block_f);
    MPI_Type_commit(&global_block_f);

    // For visualization
    int rec_len = ITERATIONS / recs;
    int nrec = 0;
    vector<double> u_global(N * N, 0);

    for (int it = 0; it < ITERATIONS; it++) {
        if (it % rec_len == 0) {
            // Fill u_global
            if (task_id != 0) {
                // Copy u to buffer
                vector<double> buffer(N * w);
                for (int i = 0; i < N; i++) 
                    for (int j = 1; j <= w; j++) 
                        buffer[i * w + (j-1)] = u[i * (w+2) + j];
                // Send buffer
                MPI_Send(buffer.data(), N*w, MPI_DOUBLE, 0, 0, 
                            MPI_COMM_WORLD);
            }
            
            // Gather in process 0
            if (task_id == 0) {
                // Copy u (of process 0) to corresponding u_global submatrix
                for (int i = 0; i < N; i++) {
                    for (int j = 1; j <= width; j++) {
                        u_global[i * N + j] = u[i * (width + 2) + j];
                    }
                }

                // Receive from tasks 1, 2, ..., (n_tasks - 2)
               for (int task = 1; task < n_tasks - 1; task++) {
                    // Recieve data in buffer
                    vector<double> buffer(N * width);
                    MPI_Recv(buffer.data(), N * width, MPI_DOUBLE, 
                                task, 0, MPI_COMM_WORLD, 
                                MPI_STATUS_IGNORE);
                    // Fill u_global with buffer data
                    int global_start_j = width * task + 1;
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < width; j++) {
                            int idx = i * N + global_start_j + j;
                            u_global[idx] = buffer[i * width + j];
                        }
                    }    
               }

               // Recieve from task (n_tasks - 1)
               vector<double> buffer(N * width_f);
               MPI_Recv(buffer.data(), N * width_f, MPI_DOUBLE, 
                        n_tasks - 1, 0, MPI_COMM_WORLD, 
                        MPI_STATUS_IGNORE);
                // Fill u_global from buffer data
               int global_start_j = width * (n_tasks - 1) + 1;
               for (int i = 0; i < N; i++) {
                    for (int j = 0; j < width_f; j++) {
                        int idx = i * N + global_start_j + j;
                        u_global[idx] = buffer[i * width_f + j];
                    }
               }

                // Write in .vtx file
                nrec++;
                write_vtx(nrec, u_global);
            }
        }

        // Exchanging ghost columns
        // [1] Sending right ghost column to the first column of the 
        //     next process.
        if (task_id % 2 == 0) {     // to avoid deadlock
            if (task_id + 1 < n_tasks) {
                MPI_Send(&u[w], 1, column, task_id + 1, 0, 
                            MPI_COMM_WORLD);
            }
            if (task_id - 1 >= 0) {
                MPI_Recv(&u[0], 1, column, task_id - 1, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            if (task_id - 1 >= 0) {
                MPI_Recv(&u[0], 1, column, task_id - 1, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (task_id + 1 < n_tasks) {
                MPI_Send(&u[w], 1, column, task_id + 1, 0,
                            MPI_COMM_WORLD);
            }
        }

        // [2] Sending left ghost column to the last column of 
        //     the previous process.
        if (task_id % 2 == 0) {     // to avoid deadlock
            if (task_id - 1 >= 0) {
                MPI_Send(&u[1], 1, column, task_id - 1, 0, 
                            MPI_COMM_WORLD);
            }
            if (task_id + 1 < n_tasks) {
                MPI_Recv(&u[w+1], 1, column, task_id + 1, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } 
        } else {
            if (task_id + 1 < n_tasks) {
                MPI_Recv(&u[w+1], 1, column, task_id + 1, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } 
            if (task_id - 1 >= 0) {
                MPI_Send(&u[1], 1, column, task_id - 1, 0, 
                            MPI_COMM_WORLD);
            }
        }

        // Fill u_new
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j <= w; j++) {
                u_new[i * (w+2) + j] = heat_diffusion(
                                            u[i * (w+2) + j],
                                            u[(i-1) * (w+2) + j],
                                            u[(i+1) * (w+2) + j],
                                            u[i * (w+2) + (j-1)],
                                            u[i * (w+2) + (j+1)]
                                            );
            }
        }
        // If necessary, reset center value to 100
        if (is_center_contained) {
            int local_j = N / 2 - start_j + 1;
            u_new[N/2 * (w+2) + local_j] = 100.0;
        }
        swap(u, u_new);
    }

    // Fill u_global
    if (task_id != 0) {
        // Copy u to buffer
        vector<double> buffer(N * w);
        for (int i = 0; i < N; i++) 
            for (int j = 1; j <= w; j++) 
                buffer[i * w + (j-1)] = u[i * (w+2) + j];
        // Send buffer
        MPI_Send(buffer.data(), N*w, MPI_DOUBLE, 0, 0, 
                    MPI_COMM_WORLD);
    }
        
    // Gather in process 0
    if (task_id == 0) {
        // Copy u (from process 0) to corresponding u_global submatrix
        for (int i = 0; i < N; i++) {
            for (int j = 1; j <= width; j++) {
                u_global[i * N + j] = u[i * (width + 2) + j];
            }
        }

        // Receive from tasks 1, 2, ..., (n_tasks - 2)
        for (int task = 1; task < n_tasks - 1; task++) {
            // Recieve data in buffer
            vector<double> buffer(N * width);
            MPI_Recv(buffer.data(), N * width, MPI_DOUBLE, 
                        task, 0, MPI_COMM_WORLD, 
                        MPI_STATUS_IGNORE);
            // Fill u_global with buffer data
            int global_start_j = width * task + 1;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < width; j++) {
                    int idx = i * N + global_start_j + j;
                    u_global[idx] = buffer[i * width + j];
                }
            }    
        }

        // Recieve from task (n_tasks - 1)
        vector<double> buffer(N * width_f);
        MPI_Recv(buffer.data(), N * width_f, MPI_DOUBLE, 
                n_tasks - 1, 0, MPI_COMM_WORLD, 
                MPI_STATUS_IGNORE);
        // Fill u_global from buffer
        int global_start_j = width * (n_tasks - 1) + 1;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < width_f; j++) {
                int idx = i * N + global_start_j + j;
                u_global[idx] = buffer[i * width_f + j];
            }
        }

        // Write in .vtx file
        nrec++;
        write_vtx(nrec, u_global);
    }
    

    // Free datatypes and finalize
    MPI_Type_free(&column);
    MPI_Type_free(&global_block);
    MPI_Type_free(&global_block_f);
    MPI_Finalize();

    return 0;   
}