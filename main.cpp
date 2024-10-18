#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>

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

vector<string> generateIdealSpectrum(int k, int n, string DNA) {
    vector<string> idealSpectrum;
    for (int i = 0; i <n-k+1; i++) {
        string oligonucleotide = DNA.substr(i, k);
        idealSpectrum.push_back(oligonucleotide);

    }
        return idealSpectrum;
    }

int main() {
    srand(static_cast<unsigned>(time(0)));
    int n,k; // n - długosć łańcucha , k - długość oligonukleotydów
    string DNA, primer;
    vector<string> idealSpectrum;
    cout << "Podaj długość łańcucha: "<<endl;
    cin >> n;
    cout << "Podaj długosć oligonukleotydów:" << endl;
    cin >> k;

    DNA =  generateDNA(n);
    idealSpectrum = generateIdealSpectrum(k, n, DNA);
    primer = idealSpectrum[0];

    cout << DNA << endl;
    cout << primer << endl;
    for (string element : idealSpectrum) {
        std::cout << element << " ";  // Wyświetlamy każdy element
    }
    return 0;
}
