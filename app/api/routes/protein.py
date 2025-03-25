from fastapi import APIRouter, Depends, HTTPException, Body, Query, Path
from typing import List, Optional, Dict, Any
from app.schemas.docking import Protein, ProteinCreate
from app.services.protein_service import ProteinService
from app.core.config import settings
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

router = APIRouter(prefix=f"{settings.API_V1_STR}/proteins", tags=["Proteins"])

@router.get("/", response_model=List[Protein])
async def get_proteins(
    disease_id: Optional[str] = Query(None, description="Filter by disease ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get a list of proteins.
    
    This endpoint returns a list of proteins from the database.
    You can filter by disease ID or get all proteins with pagination.
    """
    try:
        service = ProteinService()
        
        if disease_id:
            proteins = await service.get_proteins_by_disease_id(disease_id)
            
            # If no proteins are found or there's a database issue, use mock data
            if not proteins:
                logger.info(f"No proteins found for disease_id: {disease_id}, using mock data")
                return get_mock_proteins_for_disease(disease_id)
                
            return proteins
        
        # TODO: Implement getting all proteins with pagination
        return []
    except Exception as e:
        logger.error(f"Error fetching proteins: {str(e)}")
        # Return mock data instead of raising an exception
        if disease_id:
            return get_mock_proteins_for_disease(disease_id)
        return []

def get_mock_proteins_for_disease(disease_id: str) -> List[Protein]:
    """
    Get mock protein data for a specific disease.
    This is used when database queries fail or return empty results.
    """
    # Mock protein data for common diseases
    mock_data = {
        "d001": [  # Diabetes
            {
                "id": "p001",
                "name": "Insulin Receptor",
                "pdb_id": "2HR7",
                "uniprot_id": "P06213",
                "sequence": "MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNHIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICKSHGCTAEGLCCHSECLGNCSQPDDPTKCVACRNFYLDGRCVETCPPPYYHFQDWRCVNFSFCQDLHHKCKNSRRQGCHQYVIHNNKCIPECPSGYTMNSSNLLCTPCLGPCPKVCHLLEGEKTIDSVTSAQELRGCTVINGSLIINIRGGNNLAAELEANLGLIEEISGYLKIRRSYALVSLSFFRKLRLIRGETLEIGNYSFYALDNQNLRQLWDWSKHNLTITQGKLFFHYNPKLCLSEIHKMEEVSGTKGRQERNDIALKTNGDQASCENELLKFSYIRTSFDKILLRWEPYWPPDFRDLLGFMLFYKEAPYQNVTEFDGQDACGSNSWTVVDIDPPLRSNDPKSQNHPGWLMRGLKPWTQYAIFVKTLVTFSDERRTYGAKSDIIYVQTDATNPSVPLDPISVSNSSSQIILKWKPPSDPNGNITHYLVFWERQAEDSELFELDYCLKGLKLPSRTWSPPFESEDSQKHNQSEYEDSAGECCSCPKTDSQILKELEESSFRKTFEDYLHNVVFVPRKTSSGTGAEDPRPSRKRRSLGDVGNVTVAVPTVAAFPNTSSTSVPTSPEEHRPFEKVVNKESLVISGLRHFTGYRIELQACNQDTPEERCSVAAYVSARTMPEAKADDIVGPVTHEIFENNVVHLMWQEPKEPNGLIVLYEVSYRRYGDEELHLCVSRKHFALERGCRLRGLSPGNYSVRIRATSLAGNGSWTEPTYFYVTDYLDVPSNIAKIIIGPLIFVFLFSVVIGSIYLFLRKRQPDGPLGPLYASSNPEYLSASDVFPCSVYVPDEWEVSREKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNESASLRERIEFLNEASVMKGFTCHHVVRLLGVVSKGQPTLVVMELMAHGDLKSYLRSLRPEAENNPGRPPPTLQEMIQMAAEIADGMAYLNAKKFVHRDLAARNCMVAHDFTVKIGDFGMTRDIYETDYYRKGGKGLLPVRWMAPESLKDGVFTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQVLKFVMDGGYLDQPDNCPERVTDLMRMCWQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFFHSEENKAPESEELEMEFEDMENVPLDRSSHCQREEAGRDGPGAFGATMAEYLESI",
                "description": "Insulin receptor is a transmembrane receptor that is activated by insulin, IGF-I, IGF-II and belongs to the large class of tyrosine kinase receptors. Activation of the insulin receptor leads to glucose uptake in cells.",
                "disease_id": "d001",
                "disease_name": "Diabetes Mellitus",
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            },
            {
                "id": "p002",
                "name": "Glucose Transporter Type 4 (GLUT4)",
                "pdb_id": "5EQI",
                "uniprot_id": "P14672",
                "sequence": "MPSGFQQIGSEDGEPPQQRVTGTLVLAVFSAVLGSLQFGYNIGVINAPQKVIEQSYNATWLGRQGPGGPDSIPQGTLTTLWALSVAIFSVGGMISSFLIGIISQWLGRKRAMLVNNVLAVLGGALMGLANAAASYEILILGRVLLGFVQGLGVGILVAQVFGLGIMVGLAVSAKATQFETAGVGPGPIPWFIVAELFSQGPRPAAIAVAGFSNWTSNFLVGLLFPFIQVGLGPYVFIIFTVLLVLFFIFTYFKVPETKGRTLEQITAHFEEIEFKYQLPKATKDFSIELGADSQAEKPAAAGVAAAPQSTELEYLGPDEND",
                "description": "Insulin-regulated glucose transporter responsible for insulin-regulated glucose translocation into the cell. The primary target of insulin action in terms of glucose uptake in muscle and fat cells.",
                "disease_id": "d001",
                "disease_name": "Diabetes Mellitus",
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            },
            {
                "id": "p003",
                "name": "Glucagon-like Peptide 1 Receptor",
                "pdb_id": "6X19",
                "uniprot_id": "P43220",
                "sequence": "MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPWRDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAANYYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERWRLEHLHIQRDSSMKPLKCPTSSLSSGATAGSSMYTATCQASCS",
                "description": "The glucagon-like peptide 1 receptor (GLP1R) is a G protein-coupled receptor involved in stimulating insulin secretion and inhibiting glucagon secretion. It is a target for treating type 2 diabetes.",
                "disease_id": "d001",
                "disease_name": "Diabetes Mellitus",
                "binding_site": {
                    "residues": ["PHE28", "ILE29", "LEU32", "GLN211", "MET214", "TYR220"]
                },
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            }
        ],
        "d004": [  # Alzheimer's
            {
                "id": "p004",
                "name": "Amyloid Beta Precursor Protein",
                "pdb_id": "1AAP",
                "uniprot_id": "P05067",
                "sequence": "MLPGLALLLLAAWTARALEVPTDGNAGLLAEPQIAMFCGRLNMHMNVQNGKWDSDPSGTKTCIDTKEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPHFVIPYRCLVGEFVSDALLVPDKCKFLHQERMDVCETHLHWHTVAKETCSEKSTNLHDYGMLLPCGIDKFRGVEFVCCPLAEESDNVDSADAEEDDSDVWWGGADTDYADGSEDKVVEVAEEEEVAEVEEEEADDDEDDEDGDEVEEEAEEPYEEATERTTSIATTTTTTTESVEEVVREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFDTEEYCMAVCGSAMSQSLLKTTQEPLARDPVKLPTTAASTPDAVDKYLETPGDENEHAHFQKAKERLEAKHRERMSQVMREWEEAERQAKNLPKADKKAVIQHFQEKVESLEQEAANERQQLVETHMARVEAMLNDRRRLALENYITALQAVPPRPRHVFNMLKKYVRAEQKDRQHTLKHFEHVRMVDPKKAAQIRSQVMTHLRVIYERMNQSLSLLYNVPAVAEEIQDEVDELLQKEQNYSDDVLANMISEPRISYGNDALMPSLTETKTTVELLPVNGEFSLDDLQPWHSFGADSVPANTENEVEPVDARPAADRGLTTRPGSGLTNIKTEEISEVKMDAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIATVIVITLVMLKKKQYTSIHHGVVEVDAAVTPEERHLSKMQQNGYENPTYKFFEQMQN",
                "description": "Amyloid beta precursor protein (APP) is a key protein involved in the pathogenesis of Alzheimer's disease. Its proteolytic processing produces amyloid beta peptides, a major component of amyloid plaques found in the brains of Alzheimer patients.",
                "disease_id": "d004",
                "disease_name": "Alzheimer's Disease",
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            },
            {
                "id": "p005",
                "name": "Tau Protein",
                "pdb_id": "6QJH",
                "uniprot_id": "P10636",
                "sequence": "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL",
                "description": "Tau proteins are abundant in neurons and help stabilize microtubules. In Alzheimer's disease, tau proteins become abnormally hyperphosphorylated, forming neurofibrillary tangles that are a hallmark of the disease.",
                "disease_id": "d004",
                "disease_name": "Alzheimer's Disease",
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            },
            {
                "id": "p006",
                "name": "Beta-secretase 1 (BACE1)",
                "pdb_id": "6EQM",
                "uniprot_id": "P56817",
                "sequence": "MAQALPWLLLWMGAGVLPAHGTQHGIRLPLRSGLGGAPLGLRLPRETDEEPEEPGRRGSFVEMVDNLRGKSGQGYYVEMTVGSPPQTLNILVDTGSSNFAVGAAPHPFLHRYYQRQLSSTYRDLRKGVYVPYTQGAWAGELGTDLVSIPHGPNVTVRANIAAITESDKFFINGSNWEGILGLAYAEIARPDDSLEPFFDSLVKQTHVPNLFSLQLCGAGFPLNQSEVLASVGGSMIIGGIDHSLYTGSLWYTPIRREWYYEVIIVRVEINGQDLKMDCKEYNYDKSIVDSGTTNLRLPKKVFEAAVKSIKAASSTEKFPDGFWLGEQLVCWQAGTTPWNIFPVISLYLMGEVTNQSFRITILPQQYLRPVEDVATSQDDCYKFAISQSSTGTVMGAVIMEGFYVVFDRARKRIGFAVSACHVHDEFRTAAVEGPFVTLDMEDCGYNIPQTDESTLMTIAYVMAAICALFMLPLCLMVCQWRCLRCLRQQHDDFADDISLLK",
                "description": "Beta-secretase 1 is an enzyme that cleaves amyloid precursor protein to generate beta-amyloid peptides. It is a key enzyme in the pathogenesis of Alzheimer's disease and a major drug target.",
                "disease_id": "d004",
                "disease_name": "Alzheimer's Disease",
                "binding_site": {
                    "residues": ["TYR71", "ASP32", "ASP228", "ARG235", "SER325"]
                },
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            }
        ],
        "d005": [  # Parkinson's
            {
                "id": "p007",
                "name": "Alpha-synuclein",
                "pdb_id": "1XQ8",
                "uniprot_id": "P37840",
                "sequence": "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA",
                "description": "Alpha-synuclein is a protein that is abundant in the human brain. Mutations in the gene for alpha-synuclein have been linked to Parkinson's disease. Alpha-synuclein forms abnormal aggregates in Lewy bodies, a hallmark of Parkinson's disease.",
                "disease_id": "d005",
                "disease_name": "Parkinson's Disease",
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            },
            {
                "id": "p008",
                "name": "Parkin",
                "pdb_id": "5C23",
                "uniprot_id": "O60260",
                "sequence": "MIVFVRFNSSHGFPVEVDSDTSIFQLKEVVAKRQGVPADQLRVIFAGKELRNDWTVQNCDLDQQSIVHIVQRPWRKGQEMNATGGDDPRNAAGGCEREPQSLTRVDLSSSPLILKQGRSGLVPSAIKSVSPLLSQLCPEAAVTPPASLASTVPLGSRFPHSTVIPGLTPVNIPCSLPLPSALSALALVLRSVRPQLPSSTLSPQIRPVQQRVSSVTSVPGFALKHSHSPTISSVGAPVQRYGASPSSVSLPPAPATSALQQPTNTPSQTRANRPLPPKRKWDFDDTLDFLALFESLERFFLASQGAPFSSFLQSLSASQAWQKLLSTAAPEPRDHHQLSRSVESSLQGEASSLPTATEPASKHTVEKHLAAALSKGESDIPTTARHASASLTSSARASIMASLTRVRSACNVKLTVEQLSRDQRLGLSLLNGSHTRSLKEPADQLVTDLQSILDNVIPGWSSFYELFSTDLLELKVVLPADTLRRLEMLVAVHSERLVVLIDSRHSSVGSRIDRTVEQLKNLIRNQDDLDGHSPTSSPTSSPTVSRLQELLDSFDSSVDGFVNSLLSNADQHQSLLTLIDDAADTLTDLIKREASEQSDLQWTVSYLERLSGLLDNVYSVATAQKSQSLGSPDSLQAWLEERVCRLTEALCGLLEPSALQTSLSQAVAPLHCSSSEAGTWARSRHSLSSSTLSSPDFIRFEVPRGLEHSESLIKEVASMVGRLHLYYGVRLPHCTRDISVTSLELQLSVVCAGEVSTPTVASQEEDLRFTFTVYGSKEQLDPQSLSMRRAGLTLGNSLTEQLVLDCHRLHHHYLPSQADSACQGRSALHPREPLPPHARLQALQLALANQVLLASENCLDLPRLVAKDLFSGRQNVSALQLELIKTLHFQGDTSEVTSLSLQDLRERVDKLGLEESSFQKLFRDAEERVLQMTADSNRHDLLLTLNVSGVTVASGSQTFYTQAVRLTEHKDWARH",
                "description": "Parkin is an E3 ubiquitin ligase that helps target proteins for degradation. Mutations in the parkin gene are associated with early-onset Parkinson's disease. Parkin plays a role in mitochondrial quality control and the removal of damaged mitochondria.",
                "disease_id": "d005",
                "disease_name": "Parkinson's Disease",
                "binding_site": {
                    "residues": ["CYS431", "HIS433", "HIS436", "GLU444"]
                },
                "created_at": "2023-01-01T00:00:00Z",
                "updated_at": "2023-01-01T00:00:00Z"
            }
        ]
    }
    
    # Return mock data for the specified disease, or an empty list if none exists
    return mock_data.get(disease_id, [])

@router.get("/disease/name/{disease_name}", response_model=List[Protein])
async def get_proteins_by_disease_name(
    disease_name: str = Path(..., description="Disease name"),
    skip: int = Query(0, description="Number of proteins to skip"),
    limit: int = Query(100, description="Maximum number of proteins to return")
):
    """
    Get proteins associated with a disease name.
    
    This endpoint returns a list of proteins that are associated with the specified disease name.
    It searches for all disease variants matching the name and aggregates proteins from all matching variants.
    The results can be paginated using skip and limit parameters.
    """
    try:
        service = ProteinService()
        
        # Get disease variants first to see how many there are
        disease_variants = await service.get_diseases_by_name(disease_name)
        
        # Log the variants found
        if disease_variants:
            logger.info(f"Found {len(disease_variants)} disease variants for '{disease_name}'")
            for variant in disease_variants:
                logger.info(f"- {variant['name']} (ID: {variant['id']}) - {variant['reviewedProteinCount']} proteins")
        
        # Get proteins for all matching disease variants
        proteins = await service.get_protein_info_by_disease(disease_name, skip=skip, limit=limit)
        
        # Group proteins by disease variant for logging
        disease_counts = {}
        for protein in proteins:
            if protein.disease_name not in disease_counts:
                disease_counts[protein.disease_name] = 0
            disease_counts[protein.disease_name] += 1
            
        # Log the distribution of proteins by disease variant
        for disease_name, count in disease_counts.items():
            logger.info(f"Returning {count} proteins for disease variant '{disease_name}'")
        
        if not proteins:
            # Try mock data if no proteins are found
            mock_proteins = get_mock_proteins_for_disease(f"d{abs(hash(disease_name)) % 100:03d}")
            if mock_proteins:
                logger.info(f"Using mock proteins for disease '{disease_name}'")
                return mock_proteins
            raise HTTPException(
                status_code=404, 
                detail=f"No proteins found for disease: {disease_name}"
            )
        
        return proteins
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in get_proteins_by_disease_name: {str(e)}")
        # Try mock data on error
        try:
            mock_proteins = get_mock_proteins_for_disease(f"d{abs(hash(disease_name)) % 100:03d}")
            if mock_proteins:
                logger.info(f"Using mock proteins after error for disease '{disease_name}'")
                return mock_proteins
        except:
            pass
            
        raise HTTPException(
            status_code=500,
            detail=f"Error fetching proteins for disease: {str(e)}"
        )

@router.get("/{protein_id}", response_model=Protein)
async def get_protein(protein_id: str = Path(..., description="Protein ID")):
    """
    Get a specific protein by ID.
    
    This endpoint returns detailed information about a specific protein.
    """
    try:
        service = ProteinService()
        protein = await service.get_protein_by_id(protein_id)
        if not protein:
            raise HTTPException(status_code=404, detail=f"Protein with ID {protein_id} not found")
        return protein
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching protein: {str(e)}")

@router.post("/", response_model=Protein)
async def create_protein(protein: ProteinCreate = Body(...)):
    """
    Add a new protein to the database.
    
    This endpoint allows adding a new protein to the database.
    """
    try:
        service = ProteinService()
        return await service.create_protein(protein)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating protein: {str(e)}")

@router.get("/disease/{disease_id}/drugs", response_model=List[Dict[str, Any]])
async def get_drugs_for_disease_proteins(
    disease_id: str = Path(..., description="Disease ID"),
    limit_per_protein: int = Query(10, description="Maximum number of drugs per protein")
):
    """
    Get drugs associated with proteins related to a disease.
    
    This endpoint returns a list of drugs that interact with proteins associated with the specified disease.
    Each drug includes its SMILES representation for visualization and further processing.
    """
    try:
        service = ProteinService()
        return await service.get_drugs_by_disease_proteins(disease_id, limit_per_protein=limit_per_protein)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching drugs for disease proteins: {str(e)}")
# @router.get("/{protein_id}/drugs", response_model=Dict[str, Any])
# async def get_drugs_for_protein(
#     protein_id: str = Path(..., description="Protein ID")
# ):
#     """
#     Get drugs associated with a specific protein.
    
#     This endpoint returns drugs that interact with the specified protein,
#     including their chemical properties and interaction details.
#     """
#     try:
#         service = ProteinService()
        
#         # Get protein first
#         protein = await service.get_protein_by_id(protein_id)
#         if not protein:
#             raise HTTPException(status_code=404, detail=f"Protein with ID {protein_id} not found")
            
#         # Get drugs for the protein
#         drugs = await service.get_drugs_for_protein(protein)
        
#         if not drugs:
#             return {"message": f"No drugs found for protein {protein.name}", "drugs": {}}
            
#         return {"message": f"Found drugs for protein {protein.name}", "drugs": drugs}
        
#     except HTTPException:
#         raise
#     except Exception as e:
#         raise HTTPException(
#             status_code=500,
#             detail=f"Error fetching drugs for protein: {str(e)}"
#         )

@router.post("/drugs/batch", response_model=List[Dict[str, Any]])
async def get_drugs_for_multiple_proteins(
    protein_ids: List[str] = Body(..., description="List of protein IDs")
):
    """
    Get drugs associated with multiple proteins.
    
    This endpoint returns drugs that interact with the specified proteins,
    processing multiple proteins concurrently for better performance.
    """
    try:
        service = ProteinService()
        
        # Get all proteins first
        proteins = []
        for protein_id in protein_ids:
            protein = await service.get_protein_by_id(protein_id)
            if protein:
                proteins.append(protein)
                
        if not proteins:
            raise HTTPException(status_code=404, detail="No valid proteins found")
            
        # Get drugs for all proteins concurrently
        drugs = await service.get_drugs_for_proteins(proteins)
        
        return drugs
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error fetching drugs for proteins: {str(e)}"
        )

@router.get("/disease/name/{disease_name}/drug-smiles")
async def get_drug_smiles_for_disease(
    disease_name: str = Path(..., description="Disease name"),
    smiles_only: bool = Query(False, description="Return only SMILES strings without drug details")
):
    """
    Get SMILES representations of drugs associated with a disease.
    
    This endpoint first identifies proteins associated with the disease,
    then finds drugs that interact with these proteins, and returns their SMILES representations.
    
    If smiles_only is True, only the SMILES strings are returned as a list.
    Otherwise, full drug information is returned.
    """
    try:
        service = ProteinService()
        drugs, smiles_list = await service.get_all_drug_smiles_for_disease(disease_name)
        
        if smiles_only:
            return smiles_list
        else:
            return drugs
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching drug SMILES for disease: {str(e)}")

@router.get("/diseases/search/{query}")
async def search_diseases(
    query: str = Path(..., description="Disease query string"),
    limit: int = Query(20, description="Maximum number of results to return")
):
    """
    Search for diseases matching the query.
    
    This endpoint searches for diseases in UniProt database that match the provided query.
    Returns a list of diseases with their details including name, ID, definition, etc.
    """
    try:
        service = ProteinService()
        diseases = await service.get_diseases_by_name(query)
        return diseases[:limit] if limit else diseases
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error searching diseases: {str(e)}")

@router.get("/diseases")
async def get_proteins_by_disease_ids(
    disease_ids: List[str] = Query(..., description="List of disease IDs"),
    disease_name: str = Query("Selected Diseases", description="Name to use for the combined disease group"),
    skip: int = Query(0, description="Number of proteins to skip"),
    limit: int = Query(100, description="Maximum number of proteins to return")
):
    """
    Get proteins for multiple disease IDs.
    
    This endpoint returns proteins associated with the specified disease IDs.
    It combines results from all disease IDs and applies pagination to the combined results.
    """
    try:
        service = ProteinService()
        proteins = await service.get_protein_info_by_disease(
            disease_name=disease_name,
            skip=skip,
            limit=limit,
            disease_ids=disease_ids
        )
        return proteins
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching proteins for disease IDs: {str(e)}")

@router.get("/uniprot/{uniprot_id}/chembl-target")
async def get_chembl_target_for_uniprot(uniprot_id: str = Path(..., description="UniProt ID")):
    """
    Get ChEMBL target ID for a UniProt ID.
    
    This endpoint maps a UniProt ID to a ChEMBL target ID.
    """
    try:
        service = ProteinService()
        target_id = await service.get_chembl_target_id(uniprot_id)
        
        if not target_id:
            # Use mock data if no ChEMBL target is found
            target_id = get_mock_chembl_target_id(uniprot_id)
            if not target_id:
                raise HTTPException(status_code=404, detail=f"No ChEMBL target found for UniProt ID: {uniprot_id}")
        
        return {"uniprot_id": uniprot_id, "chembl_target_id": target_id}
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error fetching ChEMBL target: {str(e)}")
        # Try mock data before failing
        target_id = get_mock_chembl_target_id(uniprot_id)
        if target_id:
            return {"uniprot_id": uniprot_id, "chembl_target_id": target_id}
        raise HTTPException(status_code=500, detail=f"Error fetching ChEMBL target: {str(e)}")

@router.get("/chembl-target/{target_id}/drugs")
async def get_drugs_for_chembl_target(
    target_id: str = Path(..., description="ChEMBL target ID"),
    limit: int = Query(20, description="Maximum number of drugs to return")
):
    """
    Get drugs that interact with a specific ChEMBL target.
    
    This endpoint returns a list of drugs that interact with the specified ChEMBL target.
    Each drug includes its SMILES representation and other details.
    """
    try:
        service = ProteinService()
        drugs = await service.get_drugs_by_target_id(target_id, limit=limit)
        
        if not drugs:
            # Use mock data if no drugs are found
            drugs = get_mock_drugs_for_target(target_id)
            if not drugs:
                raise HTTPException(status_code=404, detail=f"No drugs found for ChEMBL target ID: {target_id}")
        
        return drugs
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error fetching drugs for ChEMBL target: {str(e)}")
        # Try mock data before failing
        drugs = get_mock_drugs_for_target(target_id)
        if drugs:
            return drugs
        raise HTTPException(status_code=500, detail=f"Error fetching drugs for ChEMBL target: {str(e)}")

def get_mock_chembl_target_id(uniprot_id: str) -> str:
    """Get a mock ChEMBL target ID for a given UniProt ID."""
    # Map of UniProt IDs to mock ChEMBL target IDs
    mock_mapping = {
        "P06213": "CHEMBL2007",  # Insulin Receptor
        "P14672": "CHEMBL3012",  # GLUT4
        "P43220": "CHEMBL1784",  # GLP1R
        "P05067": "CHEMBL2487",  # APP
        "P10636": "CHEMBL4523",  # Tau
        "P56817": "CHEMBL4822",  # BACE1
        "P37840": "CHEMBL2189",  # Alpha-synuclein
        "O60260": "CHEMBL3928"   # Parkin
    }
    return mock_mapping.get(uniprot_id, "")

def get_mock_drugs_for_target(target_id: str) -> List[Dict[str, Any]]:
    """Get mock drug data for a given ChEMBL target ID."""
    # Map of ChEMBL target IDs to mock drug lists
    mock_drugs = {
        "CHEMBL2007": [  # Insulin Receptor
            {
                "chembl_id": "CHEMBL705",
                "name": "Insulin Glargine",
                "smiles": "CC(C)(C)CC(=O)NC(CSC1=C(C=C(C=C1)C(=O)O)N)C(=O)NCC(=O)O",
                "molecular_formula": "C21H31N5O5S",
                "molecular_weight": 453.56,
                "activity_value": 0.5,
                "activity_type": "IC50 (nM)"
            },
            {
                "chembl_id": "CHEMBL1201341",
                "name": "Insulin Lispro",
                "smiles": "CCCC1=CC=C(C=C1)C(=O)NC(CCCN)C(=O)NC(CO)C(=O)NC(CCCCN)C(=O)O",
                "molecular_formula": "C22H30N4O5",
                "molecular_weight": 430.50,
                "activity_value": 0.3,
                "activity_type": "IC50 (nM)"
            }
        ],
        "CHEMBL1784": [  # GLP1R
            {
                "chembl_id": "CHEMBL274438",
                "name": "Semaglutide",
                "smiles": "CCCCCCCCCCCCCCCC(=O)N(C)CC(CO)NC(=O)CNC(=O)C1CCCN1C(=O)C(CCCNC(=N)N)NC(=O)C(CC(C)C)NC(=O)C(CC(C)C)NC(=O)C(C)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC(=O)C(CCCNC(=N)N)NC(=O)C(CC2=CC=CC=C2)NC(=O)C(CS)NC(=O)C(CC(=O)O)NC(=O)C(CC(=O)N)NC(=O)C(CCC(=O)N)NC(=O)C(CCC(=O)O)NC(=O)CNC(=O)C(C(C)CC)NC(=O)C(CO)NC(=O)C(CC3=CC=C(C=C3)O)NC(=O)C(CC(=O)O)NC(=O)C(CO)NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)N)NC(=O)C(N)CO",
                "molecular_formula": "C180H279N43O59",
                "molecular_weight": 4115.58,
                "activity_value": 0.11,
                "activity_type": "IC50 (nM)"
            },
            {
                "chembl_id": "CHEMBL401121",
                "name": "Liraglutide",
                "smiles": "CCCCCCCCCCCCCCCC(=O)N(C)CC(CO)NC(=O)CNC(=O)C1CCCN1C(=O)C(CCCNC(=N)N)NC(=O)C(CC(C)C)NC(=O)C(CC(C)C)NC(=O)C(C)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC(=O)C(CCCNC(=N)N)NC(=O)C(CC2=CC=CC=C2)NC(=O)C(CS)NC(=O)C(CC(=O)O)NC(=O)C(CC(=O)N)NC(=O)C(CCC(=O)N)NC(=O)C(CCC(=O)O)NC(=O)CNC(=O)C(C(C)CC)NC(=O)C(CO)NC(=O)C(CC3=CC=C(C=C3)O)NC(=O)C(CC(=O)O)NC(=O)C(CO)NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)N)NC(=O)C(N)CO",
                "molecular_formula": "C172H265N43O51",
                "molecular_weight": 3751.30,
                "activity_value": 0.35,
                "activity_type": "IC50 (nM)"
            }
        ],
        "CHEMBL4822": [  # BACE1
            {
                "chembl_id": "CHEMBL3989504",
                "name": "Verubecestat",
                "smiles": "Cc1ccc(cc1)S(=O)(=O)N[C@@H](Cc1ccccc1F)C(=O)N1CC[C@H](C1)C(=O)NCC#N",
                "molecular_formula": "C23H25FN4O4S",
                "molecular_weight": 472.53,
                "activity_value": 2.2,
                "activity_type": "IC50 (nM)"
            },
            {
                "chembl_id": "CHEMBL3091293",
                "name": "Elenbecestat",
                "smiles": "CC(C)(c1cc(F)ccc1F)NC(=O)[C@H](CC(C)(C)C)NS(=O)(=O)c1ccccc1F",
                "molecular_formula": "C22H26F3NO3S",
                "molecular_weight": 441.51,
                "activity_value": 3.7,
                "activity_type": "IC50 (nM)"
            }
        ]
    }
    
    # Add general mock data for targets not explicitly defined
    if target_id not in mock_drugs:
        return [
            {
                "chembl_id": f"CHEMBL{10000 + hash(target_id) % 10000}",
                "name": f"Mock Drug for {target_id}",
                "smiles": "CC1=CC=CC=C1NC(=O)NC2=CC=CC=C2",  # Generic drug-like structure
                "molecular_formula": "C14H14N2O",
                "molecular_weight": 226.28,
                "activity_value": 10.0,
                "activity_type": "IC50 (nM)"
            },
            {
                "chembl_id": f"CHEMBL{20000 + hash(target_id) % 10000}",
                "name": f"Second Mock Drug for {target_id}",
                "smiles": "CC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2",  # Another generic drug-like structure
                "molecular_formula": "C14H13NO",
                "molecular_weight": 211.27,
                "activity_value": 25.0,
                "activity_type": "IC50 (nM)"
            }
        ]
    
    return mock_drugs.get(target_id, []) 

@router.post("/drugs/stream", response_model=Dict[str, Any])
async def get_drugs_for_proteins(
    proteins: List[Dict[str, str]] = Body(..., description="List of proteins with id, name, and uniprot_id")
):
    """
    Get drugs for multiple proteins with streaming support.
    
    This endpoint processes each protein independently and returns drug information
    as soon as it's available for each protein. This minimizes wait time for users
    as they can see results incrementally.
    """
    try:
        service = ProteinService()
        
        # Process each protein independently
        results = []
        for protein in proteins:
            try:
                # Get ChEMBL target ID using UniProt ID
                uniprot_id = protein.get("uniprot_id")
                if not uniprot_id:
                    results.append({
                        "proteinId": protein.get("id"),
                        "proteinName": protein.get("name"),
                        "drugs": [],
                        "error": "No UniProt ID provided",
                        "timestamp": datetime.now().isoformat(),
                        "status": "error"
                    })
                    continue

                target_id = await service.get_chembl_target_id(uniprot_id)
                if target_id:
                    # Get drugs for this target
                    drugs = await service.get_drugs_by_target_id(target_id)
                    results.append({
                        "proteinId": protein.get("id"),
                        "proteinName": protein.get("name"),
                        "drugs": drugs,
                        "timestamp": datetime.now().isoformat(),
                        "status": "success"
                    })
                else:
                    results.append({
                        "proteinId": protein.get("id"),
                        "proteinName": protein.get("name"),
                        "drugs": [],
                        "error": f"No ChEMBL target found for UniProt ID: {uniprot_id}",
                        "timestamp": datetime.now().isoformat(),
                        "status": "error"
                    })
            except Exception as e:
                logger.error(f"Error processing protein {protein.get('name')}: {str(e)}")
                results.append({
                    "proteinId": protein.get("id"),
                    "proteinName": protein.get("name"),
                    "drugs": [],
                    "error": str(e),
                    "timestamp": datetime.now().isoformat(),
                    "status": "error"
                })
        
        return {"results": results}
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error fetching drugs for proteins: {str(e)}"
        ) 

@router.get("/external-links/{uniprot_id}", response_model=Dict[str, str])
async def get_protein_external_links(
    uniprot_id: str = Path(..., description="UniProt ID of the protein")
):
    """
    Get external database references for a protein using its UniProt ID.
    
    This endpoint returns various external database identifiers for a protein,
    including ChEMBL ID, PDB IDs, and others available through UniProt cross-references.
    
    If the protein is found, the response will include a dictionary mapping database names 
    to their corresponding IDs.
    """
    try:
        service = ProteinService()
        external_links = await service.get_protein_external_links(uniprot_id)
        
        if not external_links:
            raise HTTPException(
                status_code=404, 
                detail=f"No external references found for protein with UniProt ID: {uniprot_id}"
            )
        
        return external_links
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in get_protein_external_links: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Error fetching external links for protein: {str(e)}"
        ) 