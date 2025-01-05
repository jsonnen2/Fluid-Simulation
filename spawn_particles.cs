using UnityEngine;
using System.IO;
using System.Linq;

public class spawn_particles : MonoBehaviour
{
    [SerializeField] 
    private GameObject prefab; // Assign in Unity.
    private GameObject[] objs;
    private Material[] materials;
    [SerializeField] 
    private string filepath = "Assets/Saved_Fluid_Sims/splish_splash.txt"; // Default value. Assign in Unity.
    private StreamReader reader;
    private bool isRunning = true;
    private float next_read_time = 0f;
    [SerializeField] 
    private int num_particles = 10000; // Assign in Unity.

    void Start() {
        // Try to load file
        try {
            reader = new StreamReader(filepath);
        }
        catch (FileNotFoundException) {
            Debug.LogError($"Could not find file at {filepath}");
            isRunning = false;
            this.enabled = false;
        }
        
        // create GameObjects
        objs = new GameObject[num_particles];
        materials = new Material[num_particles];

        for (int i = 0; i < num_particles; i++) {
            objs[i] = Instantiate(prefab, Vector3.zero, Quaternion.identity);
            
            // Get and store material reference-- ClaudeAI
            Renderer renderer = objs[i].GetComponent<Renderer>();
            materials[i] = new Material(renderer.material);  // Create unique material instance
            renderer.material = materials[i];
        }

        // Spawn all .obj scene objects
        string obj_folder = Path.ChangeExtension(filepath, null);
        string[] folder_of_obj = Directory.GetFiles(obj_folder, "*.obj");
        foreach (string obj_file in folder_of_obj) {
            string obj_file_parsed = obj_file.Replace('\\', '/'); // Normalize path separators
            GameObject modelPrefab = UnityEditor.AssetDatabase.LoadAssetAtPath<GameObject>(obj_file_parsed);

            string grab_the_filename = Path.GetFileName(obj_file_parsed);
            string object_type = grab_the_filename.Substring(0, grab_the_filename.IndexOf('-'));
            if (object_type.Equals("inside_box")) continue;

            if (modelPrefab != null) {
                // Reflect the GameObject across the X-axis
                Vector3 scale = modelPrefab.transform.localScale;
                scale.x *= -1; // Ensure the scale is negative for the X-axis
                modelPrefab.transform.localScale = scale;
                Instantiate(modelPrefab, Vector3.zero, Quaternion.identity);
            }
            else {
                Debug.LogError($"Failed to load model at path: {obj_file}");
            }
        }
    }
    void Update() {
        if (!isRunning || reader == null) return;
        if (Time.time <= next_read_time) return; // Check if enough time has elasped between steps
        next_read_time = Time.time + 0.1f; // TODO: print time b/w frames

        // Read and parse 
        string line = reader.ReadLine();
        if (string.IsNullOrEmpty(line)) {
            // End the program.
            Debug.Log("End of file reached.");
            reader?.Close(); // TODO: claude says ?
            reader = null;
            isRunning = false;
            this.enabled = false;
            return;
        }

        var (position, gradient) = parse_line(line);
        
        // edit GameObjects
        for (int i = 0; i < position.Length; i++) {
            objs[i].transform.position = position[i];
            
            Color newColor = Color.Lerp(Color.red, Color.blue, gradient[i]);
            materials[i].color = newColor;
        }

        // TODO: place .obj objects from the scene
    }
    (Vector3[], float[]) parse_line(string line) {
        // Method written by ClaudeAI.
        // Remove outer brackets and split into vector strings
        line = line.Split("[[")[1].TrimEnd(']', ']');
        string[] vectorStrings = line.Split("], [");
        
        // Pre-allocate arrays
        int count = vectorStrings.Length;
        Vector3[] vectors = new Vector3[count];
        float[] wComponents = new float[count];
        
        // Parse each vector string
        for (int i = 0; i < count; i++) {
            var numbers = vectorStrings[i].Split(", ")
                                        .Select(float.Parse)
                                        .ToArray();
            
            vectors[i] = new Vector3(numbers[0], numbers[1], numbers[2]);
            wComponents[i] = numbers[3];
        }
        
        return (vectors, wComponents);
    }
    void OnDestroy()
    {
        // Clean up materials when the script is destroyed
        if (materials != null)
        {
            foreach (Material mat in materials)
            {
                if (mat != null)
                {
                    Destroy(mat);
                }
            }
        }
    }
}
