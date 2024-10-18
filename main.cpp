#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;


string generateDNA(int n) {
    string DNA;
for (int i = 0; i < n; i++) {

    const char nucleotides[] = {'A', 'C', 'T', 'G'};
    char generatedNucleotide = nucleotides[rand() % 4];
    DNA += generatedNucleotide;
}
return DNA;
}

int main() {
    int n,k; // n - długosć łańcucha , k - długość oligonukleotydów
    string DNA;
    vector<string> idealSpectrum;
    cout << "Podaj długość łańcucha: "<<endl;
    cin >> n;
    cout << "Podaj długosć oligonukleotydów:" << endl;
    cin >> k;

   DNA =  generateDNA(n);

    cout << DNA << endl;
    return 0;
}
