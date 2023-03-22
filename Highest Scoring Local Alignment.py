import itertools
import numpy as np

string1 = "LLMSQAMSGAFPGAEDWAWEFFCFENYMWFSKVLQVRCKSFTWGTWIVVQPDTGYGFRNGGGCIGVTKLAIREMSMQFQVPGCVCHCMVNGSAKESWGWPVDQSHWSHSHCDRKHAARWHSSAHFNHVNCRFVPLVMYYRWSTWMPCNTMRQQGAGCIRCPFMIPHGRKLRTARKYVYWPKQNRYDPLKGGVKGWLHFHEWNHNALMASHQVVWPHYDHLHCTHKGREQVFMLYKLVFKAVNQIMALQYVMPSHHPDWCDSKFKQTPQCYMGFWETQILVRFSTDPTMFIRHRYEYIIEMFCANTRILIESNHNIVKEFIYSGAKDKNTFCDVKTEFQFPNWQWLCISADKQSFTPEQIGNFSITEVVCDIMDWYEWSCSTNRDIVMYSCITKWKEREEQLGLGVAMFQNGHIQAAGMGPFSTSHTSGHNGHSQKRTQISVSCCDDWDKIHNMCTNVGHHLEEKTWLYCWIVPAAMLEQNMGVWKWTVMEIHSYWYERYHDYKIEYSLIFWPTCCIPTESRYNPCEHQGYFLSCYCSRVCIHWQAECQMPQTSWAYYKAEQCPKECYTTEQDVACKVCEMESRSNDIYNMSRVGLCYVDCTVWIVQHWPDPGMDFIHQANTCCGEAHLNIQWGRILKTTNMEGHYNPRLNANDHLNRPLFGYQFMHPLGQKTRMEQERGHEPDIGSDQWKSCDGHIWSMQCSIIAKFVDMHAKRCKNWHRFSAFHCCEFVSGQRFHLYGDKHTLFHEIKLFYDIHHHPMCFAAAQRDMRMFLRTERILIYQQCLGRCHCSMCMMCEETADMWVRTSMSACQAQVLGIQSIEEPYAYYMGHTWLDREAMVDGLEGTYSHPVPYMRRCDYGGGHCMYWEVHMLSVVPYVNCRWFVWTTNKSTRWPWLYTPADGAAPNHGHKHTPCPNKKVDHVNNYWLYLAVQGYFSDFGMGEMKYIWKHEWAVQHVSHNQGKYRMAKSEGGVTQQISTEVHLLPHKVFSMLRGFLALAYISLCAPCFCASAPMKVITCKQVGNMGCGYFRWASEMADSVQIYLDTPQSGTANFSDHCNMKVGGENIQVQKLFQLEVKFVNGSSDQPEGCRFDPTSGATCWKQRQKGWCAQGCSLTRAFFINTFTAQDNTSHWMHDPCIDTQFMLRVHDRIWGFIKDNPWEKSRSFIGRNIQCRFQTSPALYVITRHWQWDCFTENYGYTQGARTYPNHANFYGLSDYYGNVDTMMHPSNVYENVVMEFYISCSEDKDPKDRSMPGETRCALLACMTINHPIYLEWPARDYYDCYTREYRDWHMPHCGFFDFLQQCNDETFWRNMPNNKGRTSMVRWKVKRERGWQIVVPQCNVCIAECAQHLVAPNNIKATAGEYMKSNGAIPFPMIDGGQHEWRRKAESPHSTWPTMKCTFWTDMEPPEMKVENNTRYFYSASEQRSKLIENMWLCFQMPSCLRFQHDSIDRILCEDENEGWFFLLALWPWMMPRQLDMRCAHTMNFRMKMMSDAEASGFHAGYNACGEFTNNENSRCDRISERSVSIWRPHYSYVSFDPWSRVWAYYYQKHMDAEFFYDNTFTCIIGNIKWDNCICYVSHCFKREMLPRHQTWAIGHVCNMTWCKGFCWQVWEEPVLCSAGLYPGECPESAICPQFMNRTLSDNKDDHVYLHYMVKNAGNRPWIHCSQSPWRAHIAGMEYGEVSDFDPKHRMFSFPEKKMPNSKDIDRWYTKQRGQIQKLYACSYRIYPENSSDRPWRPFGCRIKLLRQPEGVRPGYCKVLHIKHSVMYKPGVFFHWCGYRITDAPFWICVFWHLNAMHFFMLLFHRQYKAREMVMKSRTHEADDRARTGAIDKMHKRECHSCNKKEEPKHYRVSVEPPEFVFGGPVVQCCVQANCVEVLNCQCREHCACQPWSEFPRFHTRGACNNTWLHIGSTDWFLHFECWCYNCAEVNQRPNNLRPKSSTMAVQNTSAWPVPWPQHPTVYDYACRINAQISHSCQVHIDAMLRHVAPGWFKFYMPHPEIHVPAPCDPIDMQNFETEVIYDEIRWSDIRKMLALGQSRNGDHCSAGYKPIYIQGCNISNAPGRSQCVYICHIGNLVGNVLNFANRDDHTLMGFWHFEDCCYVSSSEYHSNRNDFGLVTMNMVNYGQCRDNWQACWEQDFLFHVTIYGLYQVMEMLLCNLYQYKCHPKKWQGMRRFVCMNMRVWDTGYPMRKHPQYQFNGMGDNVDFPVGSTACEQHPQRHRRSFKISAHQNQPDHARKREECYQAKQGIFWYHQYGVDSRYNYNMCSMARMGHCPSGRMMRYMPFWRGGAGPAPLSCKDPRTSNTWGCTLRPFFRMFPFEKFRWDVDPACPLSDGSVRQRSQWGWLPRAKGGRRWQKDQANYLSLWHSTFGWMSNSDIEGPAEYVYCKYTECMLTGNECGDFADHMCIQEFVKAWNKRKETALKGWFEACWQAFFFLVIEHTMQVVMDPFKDWGQQRFCPTKTRDVTCTIQWRIDTDFHHAWIFWLDWDMKMEDKWKWQVRFYGAFEDVRAQHKASDCCMASGCNPQTAASVEEIRNHCQHQKRMSLVWVTHRIHRSTCPILHYWVGCGNPKVLVQTIPDTMCKRDGKYFGFVLHRWCHNTWPESISAMNNPKWCTNFEIIQESIGPKTGKYPADYWMEYIRFWEIECESEYIRCSAMTNAYWMNVNPYWRLQGSAQIERYKVRYGIRWQFMWLGL"
string2 = "GIQTNFDEYHVCMHIADTWPQPYERTVSPNCWIMRVDHDQCAQTTFLHIGTAQRGTKLIQFPHRPMVISEWIWGRTFVIVSNKHKKQMWTQRSRVWNNPELWNVYMLGATRQYKMMTPCAICLNRVTSKMWWWNHHVWVVFGKHAYIASFKFTLALNVATKTYVLTVAMMDRVNIGHGMMPTMSAKLAHIDCLLAHGIATQLINEVSSSFIDMLSVWNYHRQAFSTENRCQTENHMNLFWYQYQFSVDQAPVFGWNRESKRQQFSFLMDLRENYRKGCLGLRHKENDWWWGSDNLWGINNPIVDWPFPCAVDLVKYWGLVLQYAVFQWVLIWMGHWYGMMDFLDICVPARQIDWKFEFPQDLRSTPIRINELVHCDIPDGMHDKDNAGDKDARVCYVIVPRGQFPIMKCCRPWFFTFDFRRLEIKKDNMGQSQDWTVTDHAFFPYGEHLDTTVAIMGHQNHFFDEMRIVQCDIAYMHLVRKYNCTRETENLHMQFPMSWDHSYCPGVVEWRNRTMLGLGVPGAFTARSIDCRCNTQEVYRIAKDYPPLNNPWHGDIFMLYNHVPKRKDFQCLWVGCGRKMEAPLMQIVGVASDYWAVEDYKLWDRSLNQADSFCIYKPRGLKLSHTDRHFTRNDGPHLSASRAWKVALHFMGADCRQIRFNENRCCTCRKNLNYLITHRVIAKQAEINLWFYMWKKDGRWLGRYYHQSPTHETQYNQKQQTYMGDSRPECRVAPPNQLEACTYMYAGHTKQTPNCWFNQNYFCGMARNVHRQHAMIWDLCYVSVAQRANMLLREQQVCDKQFDMRYRKKYANDTMCHMVPYGRRCDYGGGHCMYWEVHMLVNCRWFVWNIKTMCCDAIWIFWLYTPADGAAPNHGHKHTPCKRLATNEDYTHNSVNNYWLYLAMIVVKQKIQWGERMIYFSDFGCKHEWAVQHVSHNLGKYRMAKSEPGVTQQISTEVHLLPHEACPYWKPGHGFSMLRALAYISLCAPCMCASQPGQFRWASEIFYLDTPQSGTDHCNMKWWILDCGGQVQVKGVMGSSNVPQVQLVNWYGCRFDWFLRTSGQKPWCAQGCSLTRAFFINTFTCRLSHHWMHDPCIDGSACGCQLQFMLRVHDRIWGFIEHCILDTGVNPWEKSRSFIRITRHWQCFTEIMDNSYSIPYGYTMWEHNWHHTGALDTIPNHANVRMEFYIKHYEPAKDPKDRSMPGETRCALLACMTIARDYYAMSHCKFICWTREYRDWHMPHCYHWFFFFDFQQYTAFFWRVMPNNKGRTRTVRWDVKRERGWQIVVTQCNVCKAECAQMKATAGEYMKSNGAIPFPMIWGDQHEWRRKAESPHSTWITMFWVGWDMNESRNNTRYFYSNMASELVWRSKLIENFWLCFQCKPSCLSFQHDIIDRILCEWIGYGSIMQENEGWFFPLAWPWMMPYQLDQACAHTMNFRNKMACGEFTNNETSRCDRFECYNHTEPQSIWRPAYSYSVDPRVWAYYYHMDAEFFYEAAICHPSNTFCICYVSHMSENTEVHCKREMLPRVTMHQTWAIGVCNMTWCKGFCWQVKAICWPPSKEEPVLCPYAGLYPGECPESKICPQFMNRTLWMNKWDHVYTKREPFRHYMVKPAGNNVGDRPWIDCSQSPWRAHIMDELNPGMDYGEVSKHRMFPEKKMPVQAELPFFDIDRWYDKHAIGAGKKQRGQPQKMGLPIQMMEYALPHAESSDRPWRPFGCRIKLLRRDDFHMPEGVFLGYCKVLHIKHSVMYHVIIAFPDIVSFLGVFFFALRWEWGEETYQMPEVITDAPFWICVFWHLNAMHFFMLGDAFHRQYKAREMQMKSRADDRARTGMIDKMGRECHSCWKFMSPGTKEEPETYRVSVEPPEFVFGGPVVQWCVQSQVANNCQPCFAREHCACQPWSEFFMIAQKRKYYTFYHVAYWLLDVVKKYLLGVATQHWVQGMRDMRNPLSEKARNFVSTTMWVIMWVYSIEHYWTLSLTFQWMEHLLANPMQMWTFFPSQMCVAQIKQMNKKHLVYVCVTHCNTPESRHICKMENMTGNHIIFGFAYATCPGDARADMSLLFVRHQWSRDPWWDQRFEHVCYWCWKGAGCIIAMSNCCDVCGYVLPEVYGNAWAVFWWQYIRIMFVAPPTPPQEQQRMVWKFWKPFHNDHMGLCEELFTWMNWCYGNNNPVHAGVMGDKNTKMFIRRLKRHEVTCGPPMEEVKFGCNITSKVRPWWQWKMECWNPKENMRMTIHHVITASQRFATTQQECWHEASPPWRMHDNEIQTDQLIEPLFIVALHNGPNRLMNWNMGISFLDTPVLAQAIYWIVQFRNLWYTYGERTHAANCFIQGMRCWEMRCMDMYSEPTKPPIPHILLPIIMMSAYAWANHDTVNPFENQTKVQPAEHTKMEHHWCMRGTNYSFKHHMFEAESVYQIIKNWTMIDMTQARHYNCLSRCHEGCKFEMWIKEAPIASATTGLDVNIIVVYDWRPSSPMKLGYKSRCQETSKFKDKKAVPDTDHGSKLAWKSWSWQPLSLISIKPGVNNPPGGHFGRCIWEEEYQTCAQKRCWIQITDKDSRCYHCYTQRIICGYTWDRTPLWKNTVCKWEEQYMQHIEFNFSGWWPASDTMMQYGIMPNAPWLMQMPGQQTQCCNRNRYFLPIMGADETGLHVFHRGWWDWFAPTYRPIGTDHRAEWGRFFNLLS"

sample_string1 = "MEANLY"
sample_string2 = "PENALTY"
penalty = 5

from Bio.Align import substitution_matrices
sub_mat = substitution_matrices.load("PAM250") 


def highest_scoring_local(string1, string2, penalty):
    alignment1 = ""
    alignment2 = ""
    dag_array = np.full((len(string2)+1, len(string1)+1), 0)
    bt_array = np.full((len(string2)+1, len(string1)+1), 0)
    
    #initialize first value
    dag_array[0,0] = 0
    
    #initialize values for first row in dag_array
    for i in range(1, len(string1)+1):
        dag_array[0, i] = max(0, dag_array[0,i-1]-penalty)
        bt_array[0, i] = 2
        
    #initialize values for first column in dag_array
    for j in range(1, len(string2)+1):
        dag_array[j, 0] = max(0, dag_array[j-1, 0]-penalty)
        bt_array[j, 0] = 1
    
    for i in range(1, len(string1)+1):
        for j in range(1, len(string2)+1):
            option1 = dag_array[j-1,i]-penalty
            
            option2 = dag_array[j,i-1]-penalty
            
            option3 = dag_array[j-1,i-1]+sub_mat[(string2[j-1], string1[i-1])]
            
            best_option = max(option1, option2, option3, 0)
            
            dag_array[j,i] = best_option
            
            if best_option == option1:
                bt_array[j, i] = 1
            
            elif best_option == option2:
                bt_array[j, i] = 2
                
            elif best_option == option3:
                bt_array[j, i] = 3
                
            else:
                bt_array[j, i] = 4
            
    print(dag_array)
    print(bt_array)
    j, i = np.unravel_index(np.argmax(dag_array, axis=None), dag_array.shape)
    orig_i = i
    orig_j = j
    while j > 0 or i > 0:  
        if bt_array[j,i] == 1:
            alignment1 = "-" + alignment1
            alignment2 = string2[j-1] + alignment2
            j -= 1
        if bt_array[j,i] == 2:
            alignment1 = string1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        if bt_array[j,i] == 3:
            alignment1 = string1[i-1] + alignment1
            alignment2 = string2[j-1] + alignment2
            j -= 1
            i -= 1
        if bt_array[j,i] == 4:
            break
        
    return (dag_array[orig_j, orig_i]), alignment1, alignment2
        

score, alignment1, alignment2 = highest_scoring_local(string1, string2, penalty)
with open("output_ba5f.txt", "w") as file:
    print(score)
    print(alignment1)
    print(alignment2)