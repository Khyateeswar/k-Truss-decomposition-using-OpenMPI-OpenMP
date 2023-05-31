#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;



// int check(vector<int> vis,int n){
//     for(int i=0;i<n;i++){
//         if(vis[i]<0){
//             return i;
//         }
//     }
//     return -1;
// }

void dfs(unordered_map<int,unordered_map<int,int>>& adj,vector<int>& vis,int i,int inde){
    if(vis[i]>0){
        return;
    }
    if(adj.find(i)==adj.end()){
        vis[i]=INT_MAX;
        return;
    }
    vis[i]=inde;
    unordered_map<int,int>::iterator it;
    for(it=(adj[i]).begin();it!=(adj[i]).end();it++){
        dfs(adj,vis,it->first,inde);
    }
    return;
}

vector<vector<int>> con(unordered_map<int,unordered_map<int,int>>& adj,int n){
    vector<vector<int>> res;
    vector<int> vis;
    int inde=1;
    for(int i=0;i<n;i++){
        vis.push_back(-1);
    }
    int t=0;
    while(t<n){
        dfs(adj,vis,t,inde);
        inde++;
        while(t<n && vis[t]>0){
            t++;
        }
    }
    // for(int i=1;i<inde;i++){
    //     vector<int> temp;
    //     for(int j=0;j<n;j++){
    //         if(vis[j]==i){
    //             temp.push_back(j);
    //         }
    //     }
    //     if(temp.size()!=0){
    //         res.push_back(temp);
    //     }
    // }
    // cout<<vis[3431]<<endl;
    // cout<<vis[3433]<<endl;
    map<int,vector<int>> m;
    for(int i=0;i<n;i++){
        if(m.find(vis[i])==m.end()){
            vector<int> temp;
            m[vis[i]]=temp;
        }
        m[vis[i]].push_back(i);
    }
    map<int,vector<int>>::iterator it;
    for(it=m.begin();it!=m.end();it++){
        if(it->first<(n+1)){
            res.push_back(it->second);
        }
    }
    return res;
}

void printv(vector<vector<int>> com){
    for(int i=0;i<com.size();i++){
        for(int j=0;j<com[i].size();j++){
            cout<<com[i][j]<<" ";
        }
        cout<<"\n";
    }
    return;
}

void printm(unordered_map<int,unordered_map<int,int>> adj){
    unordered_map<int,unordered_map<int,int>>::iterator it;
    for(it=adj.begin();it!=adj.end();it++){
        cout<<"Node "<<it->first<<endl;
        unordered_map<int,int>::iterator itr;
        for(itr = (it->second).begin();itr!=(it->second).end();itr++){
            cout<<itr->first<<" "<<itr->second<<endl;
        }
    }
    return;
}

void printa(unordered_map<int,unordered_map<int,int>> adj){
    unordered_map<int,unordered_map<int,int>>::iterator it;
    for(it=adj.begin();it!=adj.end();it++){
        cout<<"Node "<<it->first<<endl;
        unordered_map<int,int>::iterator itr;
        for(itr = (it->second).begin();itr!=(it->second).end();itr++){
            cout<<itr->first<<" ";
        }
        cout<<'\n';
    }
    return;
}

void printve(vector<int> v){
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<" ";
    }
    cout<<'\n';
}

vector<int> retadj(unordered_map<int,unordered_map<int,int>> adj,int v){
    vector<int> ans;
    unordered_map<int,int>::iterator it;
    for(it=adj[v].begin();it!=adj[v].end();it++){
        ans.push_back(it->first);
    }
    return ans;
}

void printu(unordered_map<int,int> a){
    unordered_map<int,int>::iterator it;
    for(it=a.begin();it!=a.end();it++){
        cout<<it->first<<" "<<it->second<<endl;
    }
    return;
}
void printua(unordered_map<int,int> a){
    unordered_map<int,int>::iterator it;
    for(it=a.begin();it!=a.end();it++){
        cout<<it->first<<'('<<it->second<<')'<<" ";
    }
    cout<<'\n';
    return;
}



// int get_rank(int arr[],int size,int v){
//     int sum = 0;
//     int i=0;
//     while(sum<v){

//     }
//     return i;
// }




int main(int argc, char** argv){

     MPI_Init(&argc, &argv);

    
    double start_time = MPI_Wtime();
    double end_time;

    // Get the number of processes in MPI_COMM_WORLD
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of this process in MPI_COMM_WORLD
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    string input_path ;//= "test0/test-input-0.gra";
    string output_path ;//= "test0/test-output-0.txt";
    string header_path ;//= "test0/test-header-0.dat";

    // MPI_File   file;
    // MPI_File_open(MPI_COMM_WORLD, "test5/test-output-5.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &file);
    int startk,endk;
    int p,taskid;
    //startk = 1;
    //endk = 12;

    int verbose;
    //verbose =1;
    int tot = argc;
    for(int i=1;i<tot;i++){
        string sr = argv[i];
        string f = sr.substr(0,3);
        if(f=="--i"){
            input_path = sr.substr(12,sr.size()-12);
        }else if(f=="--h"){
            header_path = sr.substr(13,sr.size()-13);
        }else if(f=="--o"){
            output_path = sr.substr(13,sr.size()-13);
        }else if(f=="--s"){
            startk = stoi(sr.substr(9,sr.size()-9));
        }else if(f=="--e"){
            endk = stoi(sr.substr(7,sr.size()-7));
        }else if(f=="--v"){
            verbose = stoi(sr.substr(10,sr.size()-10));
        }else if(f=="--p"){
            p = stoi(sr.substr(4,sr.size()-4));
        }else if(f=="--t"){
            taskid = stoi(sr.substr(9,sr.size()-9));
        }else{
            cout<<"Wrong arguments are provided at the command line"<<endl;
        }
    }

    
    
    unordered_map<int,unordered_map<int,int>> adj;
    //unordered_map<int,int> deg;

    ifstream ifs;
    ifs.open (input_path, ios::in | ios::binary );
    ifs.seekg (0, ios::end);
    int leng = ifs.tellg();
    ifs.seekg (0, ios::beg);
    char* ibuf = new char[leng];
    ifs.read (ibuf,leng);
    ifs.close();

    


    int n,m;                     // n is number if nodes and m is nuber of egdes
    n = *reinterpret_cast<int *>( ibuf );
    m = *reinterpret_cast<int *>( ibuf+4 );

    if(my_rank==0){
        cout<<n<<" "<<m<<endl;
    }

    int splitvs[world_size];
    for (int i=0; i<world_size; i++) {
        splitvs[i] = n / world_size;
        if ( i < n % world_size ) splitvs[i]++;
    }

    unordered_map<int,vector<int>> proclist;
    unordered_map<int,int> proc;
    int ij = 0;
    for(int i=0;i<world_size;i++){
        vector<int> s;
        for(int j=0;j<splitvs[i];j++ ){
            s.push_back(ij);
            proc[ij]=i;
            ij++;
        }
        proclist[i]=s;
    }
    

    ifstream ifh;
    ifh.open (header_path, ios::in | ios::binary );
    ifh.seekg (0, ios::end);
    int length = ifh.tellg();
    ifh.seekg (0, ios::beg);
    char* hbuf = new char[length];
    ifh.read (hbuf,length);
    ifh.close();

    for(int i=proclist[my_rank][0];i<proclist[my_rank][0]+proclist[my_rank].size();i++){
        unordered_map<int,int> temp;
        adj[i]=temp;
    }

    for(int i=proclist[my_rank][0];i<proclist[my_rank][0]+proclist[my_rank].size();i++){
        int sta = *reinterpret_cast<int *>( hbuf+4*i );
        int en = *reinterpret_cast<int *>( hbuf+4*i+4 );
        int nod = *reinterpret_cast<int *>( ibuf+sta );
        int deg = *reinterpret_cast<int *>( ibuf+sta+4 );
        //cout<<nod<<" "<<deg<<endl;
        
        for(int j=sta+8;j<en;j=j+4){
            adj[i][(*reinterpret_cast<int *>( ibuf+j ))]=-1;
            if(proc[(*reinterpret_cast<int *>( ibuf+j ))]==proc[i]){
                adj[(*reinterpret_cast<int *>( ibuf+j ))][i]=-1;
            }
            
            //cout<<*reinterpret_cast<int *>( ibuf+j )<<" ";
        }
        //cout<<" "<<endl;
    }
    if(my_rank==0){
            cout<<"input done"<<endl;
        }

    if (my_rank == 0 ) { 
        cout << "Time taken for input: " << end_time - start_time << endl;
    }
    // if(my_rank==1){
    //     unordered_map<int,unordered_map<int,int>>::iterator it;
    //     for(it=adj.begin();it!=adj.end();it++){
    //         cout<<it->first<<endl;
    //     }
    // }
   


    if(taskid==2){
        startk = endk;
    }




    int k=startk+2;
    string s="";
    while(k<=endk+2){


        //prefiltering 
        int sent,recv;
        int tsent,trecv;
        sent=0;
        recv=0;
        tsent=0;
        trecv=0;
        int prefil=0;
        int sum=1;
        vector<int> predil;
        // for(int i=0;i<proclist[my_rank].size();i++){
        //     if( adj[proclist[my_rank][i]].size()<k-1 ){
        //         predil.push_back(proclist[my_rank][i]);
        //     }
        // }
        unordered_map<int,unordered_map<int,int>>::iterator itt;
        for(itt=adj.begin();itt!=adj.end();itt++){
            if( adj[itt->first].size()<k-1 ){
                predil.push_back(itt->first);
            }
        }
        int size,tsize;
        tsize=0;
        while(tsize>0 || tsent>trecv){

            MPI_Status status;
            int flag = 0;
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
            while(flag){
                int v;
                MPI_Status rstatus;
                MPI_Recv(&v,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,&rstatus);
                recv++;
                MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
                adj[status.MPI_TAG].erase(v);
                if(adj[status.MPI_TAG].size()==0){
                    adj.erase(status.MPI_TAG);
                }
                if(adj[status.MPI_TAG].size()<k-1){
                    predil.push_back(status.MPI_TAG);
                }
            }



            if(predil.size()>0){
                int v = predil.back();
                predil.pop_back();
                unordered_map<int,int>::iterator it;
                for(it=adj[v].begin();it!=adj[v].end();it++){
                    // MPI_Request r;

                    MPI_Status status;
                    int flag = 0;
                    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
                    while(flag){
                        int v;
                        MPI_Status rstatus;
                        MPI_Recv(&v,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,&rstatus);
                        recv++;
                        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
                        adj[status.MPI_TAG].erase(v);
                        if(adj[status.MPI_TAG].size()==0){
                            adj.erase(status.MPI_TAG);
                        }
                        if(adj[status.MPI_TAG].size()<k-1){
                            predil.push_back(status.MPI_TAG);
                        }
                    }

                    if(proc[it->first]!=my_rank){
                        MPI_Send(&v,1,MPI_INT,proc[it->first],it->first,MPI_COMM_WORLD);
                        sent++;
                    }
                    
                }
                adj.erase(v);
            }
            size = predil.size();
            if(predil.size()==0){
                MPI_Allreduce(&sent,&tsent,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&recv,&trecv,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&size,&tsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            } 
        }



        // while(prefil==0){
        //     int chan = 0;
        //     for(int i=0;i<proclist[my_rank].size();i++){
        //         MPI_Status status;
        //         int flag = 0;
        //         MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
        //         while(flag){
        //             int v;
        //             MPI_Status rstatus;
        //             MPI_Recv(&v,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,&rstatus);
        //             MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
        //             adj[status.MPI_TAG].erase(v);
        //             chan=1;
        //         }
        //         if(adj[proclist[my_rank][i]].size()>0 && adj[proclist[my_rank][i]].size()<k-1 ){
        //             unordered_map<int,int>::iterator it;
        //             for(it=adj[proclist[my_rank][i]].begin();it!=adj[proclist[my_rank][i]].end();it++){
        //                 MPI_Send(&proclist[my_rank][i],1,MPI_INT,proc[it->first],it->first,MPI_COMM_WORLD);
        //             }
        //             adj.erase(proclist[my_rank][i]);
        //             chan=1;
        //         }
        //         MPI_Barrier(MPI_COMM_WORLD);
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     MPI_Allreduce(&chan,&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        //     if(sum==0){
        //         prefil = 1;
        //     }
        // }
        // if(my_rank==proc[0]){
        //     printua(adj[0]);
        // }
        if (my_rank == 0 ) { 
            end_time = MPI_Wtime();
            cout << "Time taken for prefiltering : " << end_time - start_time << endl;
        }
        if(k==startk+2){
        if(my_rank==0){
            cout<<"calculating support"<<endl;
        }
        // calculating support
        sent=0;
        recv = 0;
        tsent=0;
        trecv=0;
        unordered_map<int,unordered_map<int,int>>::iterator itr;
        for(itr=adj.begin();itr!=adj.end();itr++){
            // if(my_rank==1){
            //     cout<<itr->first<<endl;
            // }
                vector<int> rad = retadj(adj,itr->first);
                for(int j=0;j<world_size;j++){


                    MPI_Status status;
                    int flag = 0;
                    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
                    while(flag){
                        int v = status.MPI_TAG;
                        MPI_Status ss;
                        MPI_Request rs;
                        int sc;
                        MPI_Get_count(&status, MPI_INT, &sc);
                        int l[sc];
                        MPI_Irecv(&l,sc,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,&rs);
                        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
                        unordered_map<int,unordered_map<int,int>>::iterator it;
                        for(it=adj.begin();it!=adj.end();it++){
                            if(adj[it->first].find(v)!=adj[it->first].end()){
                                int count=0;
                                for(int i=0;i<sc;i++){
                                    if( adj[it->first].find(l[i])!=adj[it->first].end()){
                                        count++;
                                    }  
                                }
                                adj[it->first][v]=count;

                            }
                        
                        }
                        MPI_Wait(&rs,&ss);
                        recv++;
                        
                    }

                    //MPI_Request r;
                    MPI_Send(&rad[0],rad.size(),MPI_INT,j,itr->first,MPI_COMM_WORLD);
                    sent++;
                }
               //MPI_Barrier(MPI_COMM_WORLD);
        }
        int sup=0;
        while(sup==0){
            MPI_Status status;
            int flag = 0;
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
            while(flag){
                int v = status.MPI_TAG;
                MPI_Status ss;
                MPI_Request rs;
                int sc;
                MPI_Get_count(&status, MPI_INT, &sc);
                int l[sc];
                MPI_Irecv(&l,sc,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,&rs);
                MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status); 
                unordered_map<int,unordered_map<int,int>>::iterator it;
                for(it=adj.begin();it!=adj.end();it++){
                    if(adj[it->first].find(v)!=adj[it->first].end()){
                        int count=0;
                        for(int i=0;i<sc;i++){
                            if( adj[it->first].find(l[i])!=adj[it->first].end()){
                                count++;
                            }
                        }
                        adj[it->first][v]=count;
                    }
                        
                }
                MPI_Wait(&rs,&ss);
                recv++;   
            }
            MPI_Allreduce(&sent,&tsent,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(&recv,&trecv,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            if(tsent==trecv){
                sup=1;
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        //  if(my_rank==0){
        //     unordered_map<int,unordered_map<int,int>>::iterator it;
        //     for(it=adj.begin();it!=adj.end();it++){
        //         cout<<it->first<<endl;
        //     }
        // }

        }
        if(my_rank==0){
            cout<<"prefiltering done"<<endl;
        }
        if (my_rank == 0 ) { 
            end_time = MPI_Wtime();
            cout << "Time taken for support: " << end_time - start_time << endl;
        }
        //filtering edges 
        //tag 0 for asking for support deletion - also msg size =neighbour+2 source, dest
        //tag 1 for decreasing support - also msg size =2 source, dest



        sent=0;
        recv = 0;
        tsent=0;
        trecv=0;
        tsize=0;

        vector<pair<int,int>> del;
        unordered_map<int,unordered_map<int,int>>::iterator itr;
        for(itr=adj.begin();itr!=adj.end();itr++){
            unordered_map<int,int>::iterator it;
            for(it=adj[itr->first].begin();it!=adj[itr->first].end();it++){
                if(itr->first>it->first && it->second<k-2){
                    pair<int,int> temp;
                    temp.first = itr->first; //source
                    temp.second = it->first; //dest
                    del.push_back(temp);
                }
            }
        }
        while(tsize>0 || tsent>trecv){
            MPI_Status status;
            int flag = 0;
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
            while(flag){
                MPI_Request ss;
                int sc;
                MPI_Get_count(&status, MPI_INT, &sc);
                int l[sc];
                MPI_Irecv(&l,sc,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,&ss);
                recv++;
                if(status.MPI_TAG==0){
                    adj[l[sc-1]].erase(l[sc-2]);
                    for(int i=0;i<sc-2;i++){
                        if(adj[l[sc-1]].find(l[i])!=adj[l[sc-1]].end()){
                            int m[2];
                            m[0]=l[i];
                            m[1]=l[sc-2];
                            int z[2];
                            z[0]=l[sc-2];
                            z[1]=l[i];
                            MPI_Request r[4];
                            MPI_Isend(&z,2,MPI_INT,proc[z[0]],1,MPI_COMM_WORLD,&r[0]);
                            sent++;
                            MPI_Isend(&m,2,MPI_INT,proc[m[0]],1,MPI_COMM_WORLD,&r[1]);
                            sent++;
                            m[0]=l[i];
                            m[1]=l[sc-1];
                            MPI_Isend(&m,2,MPI_INT,proc[m[0]],1,MPI_COMM_WORLD,&r[2]);
                            sent++;
                            z[0]=l[sc-1];
                            z[1]=l[i];
                            MPI_Isend(&z,2,MPI_INT,proc[z[0]],1,MPI_COMM_WORLD,&r[3]);
                            sent++;
                        }
                    }
                }else{
                    if(adj[l[1]][l[0]]==k-2 ){
                        pair<int,int> t;
                        t.first = l[0]; //source > destination
                        t.second = l[1];
                        del.push_back(t);
                    }else{
                        adj[l[1]][l[0]]--;
                    }
                }
                MPI_Status sw;
                MPI_Wait(&ss,&sw);
                MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
            }
            if(del.size()>0){
                pair<int,int> temp = del.back();
                del.pop_back();
                vector<int> nei = retadj(adj,temp.first);
                nei.push_back(temp.first);
                nei.push_back(temp.second);
                MPI_Request r;
                MPI_Isend(&nei[0],nei.size(),MPI_INT,proc[temp.second],0,MPI_COMM_WORLD,&r);
                sent++;
                adj[temp.first].erase(temp.second);
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            size = del.size();
            if(del.size()==0){
                MPI_Allreduce(&sent,&tsent,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&recv,&trecv,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&size,&tsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank==0){
            cout<<"filtering edges  done"<<endl;
        }
        if (my_rank == 0 ) { 
            end_time = MPI_Wtime();
            cout << "Time taken for filtering edges: " << end_time - start_time << endl;
        }
        if(my_rank==proc[50]){
            printua(adj[50]);
        }


        

        k++;
    }



    // if(my_rank>0){
    //    MPI_Send(&s[0],s.size(),MPI_CHAR,0,77,MPI_COMM_WORLD);
        
    // }

    
    // if(my_rank==0){
    //     ofstream myfile;
    //     myfile.open (output_path);
    //     myfile<<s;
    //     //MPI_File_open(MPI_COMM_WORLD, "test5/test-output-5.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &file);
    //     for(int i=1;i<world_size;i++){
    //         MPI_Probe(i,77,MPI_COMM_WORLD,&status);
    //         int count;
    //         MPI_Get_count(&status,MPI_CHAR,&count);
    //         char buf [count];
    //         MPI_Recv(&buf,count,MPI_CHAR,i,77,MPI_COMM_WORLD,&status);
    //         for(int j=0;j<count;j++){
    //             myfile<<buf[j];
    //         }
    //         //MPI_File_write_all(file, buf,l, MPI_CHAR, &status);
    //     }
    //     myfile.close();
    //     //MPI_File_close(&file);
    // }


    
    if (my_rank == 0 ) {
        end_time = MPI_Wtime();
        cout << "Time taken: " << end_time - start_time << endl;
    }
    MPI_Finalize();
 



    return 0;

}