#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h> 

using namespace std;

void sup(unordered_map<int,unordered_map<int,int>> &adj,int i,int j){
    int a = adj[i].size();
    int b = adj[j].size();
    int count = 0;
    if(a<b){
        unordered_map<int,int>::iterator it;
        for(it=adj[i].begin();it!=adj[i].end();it++){
            if(adj[j].find(it->first)!=adj[j].end()){
                count++;
            }
        }
        adj[j][i]=count;
        adj[i][j]=count;
        return;
    }
    unordered_map<int,int>::iterator it;
    for(it=adj[j].begin();it!=adj[j].end();it++){
        if(adj[i].find(it->first)!=adj[i].end()){
            count++;
        }
    }
    adj[j][i]=count;
    adj[i][j]=count;
    return;

}

void decsup(unordered_map<int,unordered_map<int,int>> &adj,int i,int j,vector<pair<int,int>> &del,int k){
    if(adj[j].find(i)==adj[j].end()){
        return;
    }
    unordered_map<int,int>::iterator it;
    for(it=adj[j].begin();it!=adj[j].end();it++){
        if(it->first!=i && adj[i].find(it->first)!=adj[i].end()){
            //cout<<"("<<it->first<<","<<i<<")"<<"("<<adj[i][it->first]<<")"<<endl;
            //cout<<"("<<it->first<<","<<j<<")"<<"("<<adj[j][it->first]<<")"<<endl;

            adj[j][it->first]--;
            adj[it->first][j]--;
            adj[i][it->first]--;
            adj[it->first][i]--;

            if(adj[j][it->first]<k-2){
                pair<int,int> temp;
                temp.first = j;
                temp.second = it->first;
                del.push_back(temp);
            }
            if(adj[i][it->first]<k-2){
                pair<int,int> temp;
                temp.first = i;
                temp.second = it->first;
                del.push_back(temp);
            }
        }
    }
    return;
}

// correct
void prefil(unordered_map<int,unordered_map<int,int>> &adj,int n,int k,int start){// n is number of nodes and k is k-truss
    vector<int> del;
    for(int i=0;i<n;i++){
        if(adj[i].size()<k-1){
            del.push_back(i);
        }
    }
    while(del.size()>0){
        int i = del.back();
        del.pop_back();
        if(adj.find(i)!=adj.end()){
            unordered_map<int,int>::iterator it;
            for(it=adj[i].begin();it!=adj[i].end();it++){
                if(adj[it->first].find(i)!=adj[it->first].end()){
                    adj[it->first].erase(i);
                    if(adj[it->first].size()<k-1){
                        del.push_back(it->first);
                    }
                }
            }
            adj.erase(i);
        }
    }  
    for(int i=0;i<n;i++){
        if(adj[i].size()==0){
            adj.erase(i);
        }
    }       
    if(k==start){
        #pragma omp parallel for
        for(int i=0;i<n;i++){ 
            if(adj.find(i)!=adj.end() && adj[i].size()>0){
                unordered_map<int,int>::iterator it = adj[i].begin();
                for(it=adj[i].begin();it!=adj[i].end();it++){ // no change
                    if(i<it->first){
                        sup(adj,i,it->first);
                    }
                }
            }
        }
    }                                           // also calculate support
    
    return ;
}


void filedge(unordered_map<int,unordered_map<int,int>>&adj,int n,int k){// n is number of nodes and k is k-truss
    vector<pair<int,int>> del;
    for(int i=0;i<n;i++){
        //cout<<"filtering started for i = "<<i<<endl;
        if(adj.find(i)!=adj.end()){
            unordered_map<int,int>::iterator it;
            it=adj[i].begin();
            while(it!=adj[i].end()){
                if(adj[i][it->first]<k-2 && i<it->first){
                   pair<int,int> temp;
                    temp.first = i;
                    temp.second = it->first;
                    del.push_back(temp);
                    //cout<<temp.first<<","<<temp.second<<"("<<adj[temp.first][temp.second]<<")"<<endl;
                }
                it++;
            }
        }
       // cout<<"filtering done for i = "<<i<<endl;
    }
    //cout<<"first part of filtering edges done"<<endl;
    while(del.size()>0){
        pair<int,int> temp = del.back();
        del.pop_back();
        //cout<<temp.first<<","<<temp.second<<"("<<adj[temp.first][temp.second]<<")"<<endl;
        if(adj[temp.first].find(temp.second)!=adj[temp.first].end()){
            decsup(adj,temp.first,temp.second,del,k);
            adj[temp.first].erase(temp.second);
            adj[temp.second].erase(temp.first);
        } 
        
    }
    for(int i=0;i<n;i++){
        if(adj[i].size()==0){
            adj.erase(i);
        }
    }
    return ;
}

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

void printu(unordered_map<int,int> a){
    unordered_map<int,int>::iterator it;
    for(it=a.begin();it!=a.end();it++){
        cout<<it->first<<" "<<it->second<<endl;
    }
    return;
}






int main(int argc, char** argv){
    int required,provided;

     MPI_Init_thread(&argc, &argv,required,&provided);
     //omp_set_num_threads(2);

    // Get the number of processes in MPI_COMM_WORLD
    double start_time = MPI_Wtime();


    


    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of this process in MPI_COMM_WORLD
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    string input_path ;//= "test0/test-input-0.gra";
    string output_path ;//= "test0/test-output-0.txt";
    string header_path ;//= "test0/test-header-0.dat";

    // MPI_File   file;
    MPI_Status status;
    // MPI_File_open(MPI_COMM_WORLD, "test5/test-output-5.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &file);
    int startk,endk;
    //startk = 1;
    //endk = 12;

    int verbose;
    //verbose =1;
    int p,taskid;
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

    // cout<<p<<endl;
    // cout<<endk<<endl;
    // cout<<verbose<<endl;
    // cout<<taskid<<endl;

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
    

    ifstream ifh;
    ifh.open (header_path, ios::in | ios::binary );
    ifh.seekg (0, ios::end);
    int length = ifh.tellg();
    ifh.seekg (0, ios::beg);
    char* hbuf = new char[length];
    ifh.read (hbuf,length);
    ifh.close();

    for(int i=0;i<n;i++){
        unordered_map<int,int> temp;
        adj[i]=temp;
    }
    //#pragma omp parallel for
    #pragma omp parallel for 
    for(int i=0;i<n;i++){
        int sta = *reinterpret_cast<int *>( hbuf+4*i );
        int en = *reinterpret_cast<int *>( hbuf+4*i+4 );
        int nod = *reinterpret_cast<int *>( ibuf+sta );
        int deg = *reinterpret_cast<int *>( ibuf+sta+4 );
        //cout<<nod<<" "<<deg<<endl;
        for(int j=sta+8;j<en;j=j+4){
            adj[i][(*reinterpret_cast<int *>( ibuf+j ))]=-1;
            adj[(*reinterpret_cast<int *>( ibuf+j ))][i]=-1;
            //cout<<*reinterpret_cast<int *>( ibuf+j )<<" ";
        }
        //cout<<" "<<endl;
    }
    //cout<<"input done"<<endl;
    double end_time = MPI_Wtime();
    if (my_rank == 0 ) { 
        cout << "Time taken for input: " << end_time - start_time << endl;
    }

    unordered_map<int,unordered_map<int,int>> old_adj = adj;
    int workloads[world_size];
    int my_start,my_end;
    


    if(taskid==2){
        if(my_rank==0){
             my_start = endk;
            my_end = endk+1;
        }else{
             my_start = endk;
            my_end = endk;
       
        }
        
    }else{
        for (int i=0; i<world_size; i++) {
            workloads[i] = (endk-startk+1) / world_size;
            if ( i < (endk-startk+1) % world_size ) workloads[i]++;
        }
        my_start = startk;
        for (int i=0; i<my_rank; i++) {
            my_start +=	workloads[i];
        }
        my_end = my_start	+ workloads[my_rank];
    }

    int k=my_start+2;
    string s="";
    while(k<my_end+2){
   
        prefil(adj,n,k,my_start+2);
        end_time = MPI_Wtime();
        if (my_rank == 0 ) { 
            cout << "Time taken for prefiltering and support: " << end_time - start_time << endl;
        }
        filedge(adj,n,k);
        end_time = MPI_Wtime();
        if (my_rank == 0 ) { 
            cout << "Time taken for filtering edges: " << end_time - start_time << endl;
        }
        
  
        // printu(adj[3433]);
        // cout<<" "<<endl;
        // printu(adj[3431]);
        // cout<<" "<<endl;
        vector<vector<int>> com = con(adj,n);
        // cout<<"for k = "<<k-2<<endl;
        // cout<<com.size()<<endl;
        // printv(com);

        if(taskid==1){
            

            if(verbose==1){
                if(com.size()>0){
                    s=s+"1\n";
                }else{
                    s=s+"0\n";
                }
                if(com.size()>0){
                    s=s+to_string(com.size())+"\n";
                    for(int i=0;i<com.size();i++){
                        //sort(com[i].begin(),com[i].end());
                        for(int j=0;j<com[i].size();j++){
                            s=s+to_string(com[i][j])+" ";
                        }
                        s=s+"\n";
                    }
                }
            
            }else{
                if(com.size()>0){
                    s=s+"1 ";
                }else{
                    s=s+"0 ";
                }
            }
        }else{
            unordered_map<int,int> tr;
            unordered_map<int,string> str;
            #pragma parallel for
            for(int i=0;i<com.size();i++){
                string d = "";
                for(int j=0;j<com[i].size();j++){
                    tr[com[i][j]]=i;
                    d=d+to_string(com[i][j])+" ";
                }
                str[i]=d;
            }
            
            int count = 0;
            #pragma omp parallel for 
            for(int i=0;i<n;i++){
                int f=-1;
                if(tr.find(i)!=tr.end()){
                    f=tr[i];
                }
                unordered_set<int> ego;
                if(f!=-1){
                    ego.insert(f);
                }
                unordered_map<int,int>::iterator it;
                for(it=old_adj[i].begin();it!=old_adj[i].end();it++){
                    if(tr.find(it->first)!=tr.end()){
                        if(ego.find(tr[it->first])==ego.end()){
                            ego.insert(tr[it->first]);
                        }
                    }
                }
                string d = "";
                if(ego.size()>=p){
                    #pragma omp atomic update
                    count++;
                    
                    if(verbose==0){
                        d=d+to_string(i)+" ";
                    }else{
                        d=d+to_string(i);
                        d=d+'\n';
                    }
                    if(verbose==1){
                        unordered_set<int>::iterator it;
                        for(it=ego.begin();it!=ego.end();it++){
                            d=d+str[*it];
                        }
                        d=d+'\n';
                    }
                    #pragma omp critical
                    {
                        s=s+d;
                    }

                }
            }


            s=to_string(count)+'\n'+s;
        }
        end_time = MPI_Wtime();
        if (my_rank == 0 ) { 
            cout << "Time taken for preparing output: " << end_time - start_time << endl;
        }
        
        k++;
    }
    // MPI_File_write_all(file, (char*)&s[0], s.size(), MPI_CHAR, &status);
    // MPI_File_close(&file);

    if(taskid==1){
        if(my_rank>0){
            MPI_Send(&s[0],s.size(),MPI_CHAR,0,77,MPI_COMM_WORLD);
        
        }

    
        if(my_rank==0){
            ofstream myfile;
            myfile.open (output_path);
            myfile<<s;
            //MPI_File_open(MPI_COMM_WORLD, "test5/test-output-5.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &file);
            for(int i=1;i<world_size;i++){
                MPI_Probe(i,77,MPI_COMM_WORLD,&status);
                int count;
                MPI_Get_count(&status,MPI_CHAR,&count);
                char buf [count];
                MPI_Recv(&buf,count,MPI_CHAR,i,77,MPI_COMM_WORLD,&status);
                for(int j=0;j<count;j++){
                    myfile<<buf[j];
                }
                //MPI_File_write_all(file, buf,l, MPI_CHAR, &status);
            }
            myfile.close();
        //MPI_File_close(&file);
        }
    }else{
        if(my_rank==0){
            ofstream myfile;
            myfile.open (output_path);
            myfile<<s;
            myfile.close();
        }
    }
    


    end_time = MPI_Wtime();
    if (my_rank == 0 ) {
        cout << "Total Time taken: " << end_time - start_time << endl;
    }
    MPI_Finalize();
 



    return 0;

}
