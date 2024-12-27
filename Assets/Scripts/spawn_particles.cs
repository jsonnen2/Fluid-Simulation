using UnityEngine;
using System.IO;
using System.Linq;

public class spawn_particles : MonoBehaviour
{
    [SerializeField] 
    private GameObject prefab;
    [SerializeField]
    private string filepath = "Assets/Saved_Fluid_Sims/trial_1.txt";
    private StreamReader reader;
    bool isRunning = true;
    private float next_read_time = 0f;

    void Start() {
        try {
            reader = new StreamReader(filepath);
        }
        catch (FileNotFoundException) {
            Debug.LogError($"Could not find file at {filepath}");
            isRunning = false;
            this.enabled = false;
        }
    }

    void Update() {
        if (!isRunning) return;
        if (Time.time <= next_read_time) return; // Check if enough time has elasped between steps
        next_read_time = Time.time + 1f;

        // Destroy all previous particle GameObjects
        GameObject[] obj_collection = GameObject.FindGameObjectsWithTag("Particle");
        foreach (GameObject obj in obj_collection) {
            Destroy(obj);
        }

        // Read and parse 
        string line = reader.ReadLine();
        if (line != null) {
            Debug.Log(line.Length);

            // Split into individual vector strings and parse. ClaudeAI wrote this line.
            line = line.Split("[[")[1].TrimEnd(']', ']');
            Vector3[] vectors = line.Split("], [")
                .Select(vectorStr => 
                {
                    var numbers = vectorStr.Split(", ")
                        .Select(float.Parse)
                        .ToArray();
                    return new Vector3(numbers[0], numbers[1], numbers[2]);
                })
                .ToArray();
            foreach (Vector3 position in vectors) {
                // Debug.Log($"Vector3({vec.x}, {vec.y}, {vec.z})");
                GameObject obj = Instantiate(prefab, position, Quaternion.identity);
                obj.tag = "Particle";
            }
        }

        // End the program.
        else {
            Debug.Log("End of file reached.");
            reader.Close();
            reader = null;
            isRunning = false;
            this.enabled = false;
        }
    }
}
